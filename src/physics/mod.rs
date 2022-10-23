mod stream;

use crate::config::CONFIG;
use cgmath::{prelude::*, Vector2};
use rayon::prelude::{IntoParallelRefMutIterator, ParallelIterator};
use std::{iter::Iterator, mem};

const GRID_COLS: usize = CONFIG.physics.grid_cols;
const GRID_ROWS: usize = CONFIG.physics.grid_rows;

const RELAXATION_FACTOR: f32 = {
    let relaxation_time = 0.5
        + (3 * CONFIG.physics.grid_rows) as f32 / CONFIG.physics.reynolds_number_undisturbed_flow;
    1.0 / relaxation_time
};
const CELL_DEFAULT_VELOCITY: Vector2<f32> = {
    let sqrt3 = 1.732;
    Vector2::new(
        CONFIG.physics.mach_number_default_velocity.x / sqrt3,
        CONFIG.physics.mach_number_default_velocity.y / sqrt3,
    )
};

const LATTICE: [Vector2<f32>; 9] = [
    Vector2::new(-1.0, -1.0),
    Vector2::new(-1.0, 0.0),
    Vector2::new(-1.0, 1.0),
    Vector2::new(0.0, -1.0),
    Vector2::new(0.0, 0.0),
    Vector2::new(0.0, 1.0),
    Vector2::new(1.0, -1.0),
    Vector2::new(1.0, 0.0),
    Vector2::new(1.0, 1.0),
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
    ticks_elapsed: usize,
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
    pub fn total_cell_mass(&self) -> f32 {
        self.lattice.iter().flatten().flatten().sum()
    }
    pub fn ticks_elapsed(&self) -> usize {
        self.ticks_elapsed
    }
    pub fn new() -> Self {
        let velocity_distribution = equilibrium(CELL_DEFAULT_VELOCITY);
        let lattice = Box::new([[velocity_distribution; GRID_COLS]; GRID_ROWS]);

        let max_par = rayon::max_num_threads();
        let leaf_row_count = usize::max(10, (GRID_ROWS + max_par - 1) / max_par);
        let max_leaf_count = (GRID_ROWS + leaf_row_count - 1) / leaf_row_count;

        Self {
            lattice: lattice.clone(),
            blocked: Box::new([[false; GRID_COLS]; GRID_ROWS]),
            swapped_out_lattice: lattice,
            ticks_elapsed: 0,
            stream_fork_join_tree_leaf_row_count: leaf_row_count,
            stream_fork_join_scratch_rows: (0..max_leaf_count.next_power_of_two() * 2)
                .map(|_| [[0.0; 9]; GRID_COLS])
                .collect(),
        }
    }
    pub fn lbm_collide(&mut self) {
        self.lattice
            .par_iter_mut()
            .for_each(|lattice_row| lattice_row.iter_mut().for_each(collide_cell))
    }
    pub fn lbm_stream(&mut self) {
        mem::swap(&mut self.lattice, &mut self.swapped_out_lattice);
        stream::stream_row_range(
            &self.blocked,
            &*self.swapped_out_lattice,
            &mut *self.lattice,
            &mut self.stream_fork_join_scratch_rows,
            0,
            self.stream_fork_join_tree_leaf_row_count,
        );
        self.lattice.iter_mut().for_each(|row| {
            set_vel(&mut row[0], CELL_DEFAULT_VELOCITY);
            set_vel(&mut row[GRID_COLS - 1], CELL_DEFAULT_VELOCITY);
        });
        self.ticks_elapsed += 1;
    }
    pub fn get_velocity_vorticity_speed(
        &self,
        velocity: &mut [[Vector2<f32>; GRID_COLS]; GRID_ROWS],
        vorticity: &mut [[f32; GRID_COLS]; GRID_ROWS],
        speed: &mut [[f32; GRID_COLS]; GRID_ROWS],
    ) {
        for (r, c) in cells() {
            let vel = &mut velocity[r][c];
            *vel =
                momentum_at(&self.blocked, &self.lattice, (r, c)) / self.lattice[r][c].iter().sum();
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
        pub fn momentum_at(
            blocked: &[[bool; GRID_COLS]; GRID_ROWS],
            lattice: &[[[f32; 9]; GRID_COLS]; GRID_ROWS],
            (r, c): (usize, usize),
        ) -> Vector2<f32> {
            if blocked[r][c] {
                return Vector2::zero();
            }
            lattice[r][c]
                .iter()
                .enumerate()
                .fold(Vector2::zero(), |acc, (i, &val)| acc + LATTICE[i] * val)
        }
    }
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

fn equilibrium(vel: Vector2<f32>) -> [f32; 9] {
    let mut ret = [0.0; 9];
    for i in 0..9 {
        let along = LATTICE[i].dot(vel);
        ret[i] = WEIGHT[i] * (1.0 + 3.0 * along + 4.5 * along * along - 1.5 * vel.magnitude2());
        if ret[i].is_sign_negative() {
            panic!(
                "Rescaled speed {}, velocity {:?} is too high for D2Q9",
                vel.magnitude(),
                vel
            );
        }
    }
    ret
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
        .zip(LATTICE.iter())
        .map(|(&vel_mass, &velocity)| vel_mass * velocity)
        .sum();
    debug_assert!(momentum.is_finite());

    let vel = momentum / mass;
    debug_assert!(vel.is_finite());

    cell.iter_mut()
        .zip(equilibrium(vel).map(|eq| mass * eq))
        .for_each(|(f, f_eq)| *f += RELAXATION_FACTOR * (f_eq - *f));
}
