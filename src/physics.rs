use crate::config::CONFIG;
use cgmath::{prelude::*, Vector2};
use rayon::prelude::{IntoParallelIterator, IntoParallelRefMutIterator, ParallelIterator};
use std::{convert::TryInto, iter::Iterator, mem};

const SOUND_SPEED_SI: f32 = CONFIG.physics.sound_speed_si;
const KIN_VISCOSITY_SI: f32 = CONFIG.physics.kin_viscosity_si;
const DENSITY_SI: f32 = CONFIG.physics.density_si;
const GRID_COLS: usize = CONFIG.physics.grid_cols;
const GRID_ROWS: usize = CONFIG.physics.grid_rows;
const CELL_SIZE_SI: f32 = CONFIG.physics.cell_size_si;
const DEFAULT_VELOCITY_SI: Vector2<f32> = CONFIG.physics.default_velocity_si;

// Derived
const SQRT3: f32 = 1.732; // sqrt() is not a const fn
const TIME_STEP_SI: f32 = CELL_SIZE_SI / (SQRT3 * SOUND_SPEED_SI);
const RELAXATION_FACTOR: f32 =
    1.0 / (0.5 + (SQRT3 * KIN_VISCOSITY_SI) / (SOUND_SPEED_SI * CELL_SIZE_SI));
const CELLS: usize = GRID_COLS * GRID_ROWS;
const DEFAULT_VELOCITY: Vector2<f32> = Vector2::new(
    DEFAULT_VELOCITY_SI.x / (SQRT3 * SOUND_SPEED_SI),
    DEFAULT_VELOCITY_SI.y / (SQRT3 * SOUND_SPEED_SI),
);

const LATTICE: [(isize, isize); 9] = [
    (-1, -1),
    (-1, 0),
    (-1, 1),
    (0, -1),
    (0, 0),
    (0, 1),
    (1, -1),
    (1, 0),
    (1, 1),
];
const WEIGHT: [f32; 9] = {
    let a = 1.0 / 9.0;
    let b = 1.0 / 36.0;
    [b, a, b, a, 4.0 / 9.0, a, b, a, b]
};

pub struct Physics {
    lattice: Box<[[[f32; 9]; GRID_COLS]; GRID_ROWS]>,
    blocked: Box<[[bool; GRID_COLS]; GRID_ROWS]>,
    swapped_out_lattice: Box<[[[f32; 9]; GRID_COLS]; GRID_ROWS]>,
    stream_fork_join_tree_leaf_row_count: usize,
    stream_fork_join_scratch_rows: Vec<[[f32; 9]; GRID_COLS]>,
    // Measured in abstract time units. A cell has width 1 space unit and the
    // speed of sound is roughly 1 space unit per time unit.
    time_elapsed: usize,
}

impl Physics {
    pub fn blocked(&self) -> &[[bool; GRID_COLS]; GRID_ROWS] {
        &self.blocked
    }
    pub fn update_cell_blocked_status(&mut self, (r, c): (usize, usize), make_blocked: bool) {
        if self.blocked[r][c] == make_blocked {
            return;
        }
        if self.blocked[r][c] {
            self.lattice[r][c][4] = 1.0;
            set_vel(&mut self.lattice[r][c], Vector2::zero());
        }
        self.blocked[r][c] ^= true;
    }
    pub fn total_mass(&self) -> f32 {
        let mut mass = 0.0;
        for (r, c) in cells() {
            for v in 0..9 {
                mass += self.lattice[r][c][v];
            }
        }
        mass * (DENSITY_SI * CELL_SIZE_SI.powi(3))
    }
    pub fn mass_at(&self, (r, c): (usize, usize)) -> f32 {
        self.lattice[r][c].iter().sum()
    }
    pub fn momentum_at(&self, (r, c): (usize, usize)) -> Vector2<f32> {
        if self.blocked[r][c] {
            return Vector2::zero();
        }
        self.lattice[r][c]
            .iter()
            .enumerate()
            .fold(Vector2::zero(), |acc, (i, &val)| {
                acc + vec2_from_pair(LATTICE[i]) * val
            })
    }
    pub fn time_elapsed(&self) -> f32 {
        self.time_elapsed as f32 * TIME_STEP_SI
    }
    pub fn new() -> Self {
        let mut lattice = Self::new_lattice();
        let equiv = {
            let mut equiv = [0.0; 9];
            equiv
                .iter_mut()
                .zip(equilibrium(DEFAULT_VELOCITY))
                .for_each(|(arr, eq)| *arr = eq);
            equiv
        };
        lattice
            .iter_mut()
            .flat_map(|lattice_row| lattice_row.iter_mut())
            .for_each(|cell| *cell = equiv);
        let blocked = vec![false; CELLS];
        let blocked: Box<[bool; CELLS]> = blocked.into_boxed_slice().try_into().unwrap();

        let max_par = rayon::max_num_threads();
        let leaf_row_count = usize::max(10, (GRID_ROWS + max_par - 1) / max_par);
        let max_leaf_count = (GRID_ROWS + leaf_row_count - 1) / leaf_row_count;

        unsafe {
            Self {
                lattice: lattice.clone(),
                blocked: std::mem::transmute(blocked),
                swapped_out_lattice: lattice,
                time_elapsed: 0,
                stream_fork_join_tree_leaf_row_count: leaf_row_count,
                stream_fork_join_scratch_rows: (0..max_leaf_count.next_power_of_two() * 2)
                    .map(|_| [[0.0; 9]; GRID_COLS])
                    .collect(),
            }
        }
    }
    fn new_lattice() -> Box<[[[f32; 9]; GRID_COLS]; GRID_ROWS]> {
        let lattice = vec![0.0; 9 * CELLS];
        let lattice: Box<[f32; 9 * CELLS]> = lattice.into_boxed_slice().try_into().unwrap();
        unsafe { std::mem::transmute(lattice) }
    }
    pub fn lbm_collide(&mut self) {
        self.lattice
            .par_iter_mut()
            .for_each(|lattice_row| lattice_row.iter_mut().for_each(collide_cell))
    }
    pub fn lbm_stream(&mut self) {
        mem::swap(&mut self.lattice, &mut self.swapped_out_lattice);
        stream_row_range(
            &self.blocked,
            &*self.swapped_out_lattice,
            &mut *self.lattice,
            &mut self.stream_fork_join_scratch_rows,
            0,
            self.stream_fork_join_tree_leaf_row_count,
        );
        self.lattice.iter_mut().for_each(|row| {
            set_vel(&mut row[0], DEFAULT_VELOCITY);
            set_vel(&mut row[GRID_COLS - 1], DEFAULT_VELOCITY);
        });
        self.time_elapsed += 1;
    }
    pub fn get_velocity_vorticity_speed(
        &self,
        velocity: &mut [[Vector2<f32>; GRID_COLS]; GRID_ROWS],
        vorticity: &mut [[f32; GRID_COLS]; GRID_ROWS],
        speed: &mut [[f32; GRID_COLS]; GRID_ROWS],
    ) {
        for (r, c) in cells() {
            let vel = &mut velocity[r][c];
            *vel = self.momentum_at((r, c)) / self.mass_at((r, c));
            if !vel.is_finite() {
                *vel = Vector2::zero();
            }
            speed[r][c] = vel.magnitude();
        }
        for (r, c) in cells() {
            if self.blocked[r][c] {
                continue;
            }
            let mut res = 0.0;
            if r != 0 && r != GRID_ROWS - 1 {
                res -= velocity[r - 1][c].y;
                res += velocity[r + 1][c].y;
            }
            if c != 0 && c != GRID_COLS - 1 {
                res += velocity[r][c - 1].x;
                res -= velocity[r][c + 1].x;
            }
            vorticity[r][c] = if res.is_finite() { res } else { 0.0 };
        }
    }
}

const fn vec2_from_pair((x, y): (isize, isize)) -> Vector2<f32> {
    Vector2::new(x as f32, y as f32)
}

fn set_vel(cell: &mut [f32; 9], vel: Vector2<f32>) {
    let mass: f32 = cell.iter().sum();
    cell.iter_mut()
        .zip(equilibrium(vel))
        .for_each(|(cell, eq)| *cell = mass * eq);
}

pub fn cells() -> impl Iterator<Item = (usize, usize)> {
    (0..GRID_ROWS).flat_map(|r| (0..GRID_COLS).map(move |c| (r, c)))
}

fn equilibrium(vel: Vector2<f32>) -> impl Iterator<Item = f32> {
    (0..9)
        .map(move |v| {
            let along = vec2_from_pair(LATTICE[v]).dot(vel);

            WEIGHT[v] * (1.0 + 3.0 * along + 4.5 * along * along - 1.5 * vel.magnitude2())
        })
        .filter(move |f_eq| {
            if f_eq.is_sign_negative() {
                panic!(
                    "Rescaled speed {}, velocity {:?} is too high for D2Q9",
                    vel.magnitude(),
                    vel
                );
            }
            true
        })
}

fn collide_cell(cell: &mut [f32; 9]) {
    let mass: f32 = cell.iter().sum();
    debug_assert!(mass.is_finite());
    if mass.is_zero() {
        return;
    }
    debug_assert!(mass.is_sign_positive());

    let momentum: Vector2<f32> = cell
        .iter()
        .enumerate()
        .map(|(i, &p)| vec2_from_pair(LATTICE[i]) * p)
        .sum();
    debug_assert!(momentum.is_finite());

    let vel = momentum / mass;
    debug_assert!(vel.is_finite());

    cell.iter_mut()
        .zip(equilibrium(vel).map(|eq| mass * eq))
        .for_each(|(f, f_eq)| *f += RELAXATION_FACTOR * (f_eq - *f));
}
fn first_last_mut<T>(slice: &mut [T]) -> (&mut T, &mut T) {
    let (a, b) = slice.split_first_mut().unwrap();
    (a, b.last_mut().unwrap())
}
fn double_index_mut<T>(slice: &mut [T], i: usize) -> (&mut T, &mut T) {
    assert!(slice.len() > i + 1);
    let (a, b) = slice.split_at_mut(i + 1);
    (&mut a[i], &mut b[0])
}
fn triple_index_mut<T>(slice: &mut [T], i: usize) -> (&mut T, &mut T, &mut T) {
    assert!(slice.len() > i + 2);
    let (a, rest) = slice.split_at_mut(i + 1);
    let (b, c) = rest.split_first_mut().unwrap();
    (&mut a[i], b, &mut c[0])
}

fn stream_row_range(
    blocked: &[[bool; GRID_COLS]; GRID_ROWS],
    old_lattice_rows: &[[[f32; 9]; GRID_COLS]],
    new_lattice_rows: &mut [[[f32; 9]; GRID_COLS]],
    scratch_rows: &mut [[[f32; 9]; GRID_COLS]],
    base_row: usize,
    leaf_row_count: usize,
) {
    assert_eq!(old_lattice_rows.len(), new_lattice_rows.len());
    if old_lattice_rows.len() <= leaf_row_count {
        let (row_before, row_after) = first_last_mut(scratch_rows);
        row_before.fill([0.0; 9]);
        row_after.fill([0.0; 9]);
        new_lattice_rows.fill([[0.0; 9]; GRID_COLS]);
        if old_lattice_rows.len() == 1 {
            stream_row(
                blocked,
                old_lattice_rows.first().unwrap(),
                [row_before, &mut new_lattice_rows[0], row_after],
                base_row,
            );
        } else {
            let (first, second) = double_index_mut(new_lattice_rows, 0);
            stream_row(
                blocked,
                &old_lattice_rows[0],
                [row_before, first, second],
                base_row,
            );
            for i in 0..(old_lattice_rows.len() - 2) {
                let (a, b, c) = triple_index_mut(new_lattice_rows, i);
                stream_row(
                    blocked,
                    &old_lattice_rows[i + 1],
                    [a, b, c],
                    base_row + i + 1,
                );
            }
            let (secondlast, last) = double_index_mut(new_lattice_rows, new_lattice_rows.len() - 2);
            stream_row(
                blocked,
                old_lattice_rows.last().unwrap(),
                [secondlast, last, row_after],
                base_row + old_lattice_rows.len() - 1,
            );
        }
    } else {
        let mid = old_lattice_rows.len() / 2;
        let (old_left, old_right) = old_lattice_rows.split_at(mid);
        let (new_left, new_right) = new_lattice_rows.split_at_mut(mid);
        let (scratch_left, scratch_right) = scratch_rows.split_at_mut(scratch_rows.len() / 2);
        rayon::join(
            || {
                stream_row_range(
                    blocked,
                    old_left,
                    new_left,
                    scratch_left,
                    base_row,
                    leaf_row_count,
                )
            },
            || {
                stream_row_range(
                    blocked,
                    old_right,
                    new_right,
                    scratch_right,
                    base_row + mid,
                    leaf_row_count,
                )
            },
        );
        for c in 0..GRID_COLS {
            for vel in 0..9 {
                new_left.last_mut().unwrap()[c][vel] += scratch_right.first().unwrap()[c][vel];
                new_right.first_mut().unwrap()[c][vel] += scratch_left.last().unwrap()[c][vel];
            }
        }
    }
}
fn stream_row(
    blocked: &[[bool; GRID_COLS]; GRID_ROWS],
    old_row: &[[f32; 9]; GRID_COLS],
    new_rows: [&mut [[f32; 9]; GRID_COLS]; 3],
    r0: usize,
) {
    for (c0, old_cell) in old_row.iter().enumerate() {
        if blocked[r0][c0] {
            continue;
        }
        for (dr0, dc0) in LATTICE {
            let (r, c, dr, dc) = step_cell_coordinates(blocked, (r0, c0), (dr0, dc0));
            debug_assert!(r < GRID_ROWS);
            debug_assert!(c < GRID_COLS);
            debug_assert!(!blocked[r][c]);
            let old_vel = (3 * dr0 + dc0 + 4) as usize;
            let new_vel = (3 * dr + dc + 4) as usize;
            new_rows[r.wrapping_sub(r0).wrapping_add(1)][c][new_vel] += old_cell[old_vel];
        }
    }
}

fn step_cell_coordinates(
    blocked: &[[bool; GRID_COLS]; GRID_ROWS],
    (r, c): (usize, usize),
    (dr, dc): (isize, isize),
) -> (usize, usize, isize, isize) {
    // Bounce at upper and lower wall
    let (nr, dr) = {
        let nr = r as isize + dr;
        if nr == -1 {
            (0, 1)
        } else if nr == GRID_ROWS as isize {
            (GRID_ROWS - 1, -1)
        } else {
            (nr as usize, dr)
        }
    };
    // Periodic at left and right wall
    let nc = {
        let nc = c as isize + dc;
        if nc == -1 {
            GRID_COLS - 1
        } else if nc == GRID_COLS as isize {
            0
        } else {
            nc as usize
        }
    };
    // Bounce against blocked cells
    let blocked_both = blocked[nr][nc];
    let blocked_r = blocked[nr][c];
    let blocked_c = blocked[r][nc];
    if (blocked_r && blocked_c) || (blocked_both && !blocked_r && !blocked_c) {
        (r, c, -dr, -dc)
    } else if !blocked_both {
        (nr, nc, dr, dc)
    } else if blocked_c {
        (nr, c, dr, -dc)
    } else {
        (r, nc, -dr, dc)
    }
}
