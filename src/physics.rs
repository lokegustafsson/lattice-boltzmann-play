use cgmath::{prelude::*, Vector2};
use std::convert::TryInto;
use std::iter::Iterator;

// LBM restricts us to relatively low reynolds numbers. Also, since this implementation is very
// slow, we need to keep a significant fraction of the speed of sound to see any action. Hence,
// this simulation simulates a few grams of honey flowing through a tube the size of a water
// bottle at 100 m/s, in extreme slow motion

// Physical (honey, roughly)
const SOUND_SPEED_SI: f32 = 1400.0; // m/s
const KIN_VISCOSITY_SI: f32 = 0.005; // m^2/s
const DENSITY_SI: f32 = 1000.0; // kg/m^3

// Constants
pub const PHYS_COLS: usize = 200;
pub const PHYS_ROWS: usize = 120;
const CELL_SIZE_SI: f32 = 0.001; // 1 mm
const DEFAULT_VELOCITY_SI: Vector2<f32> = Vector2::new(0.0, 100.0); // 100 m/s

// Derived
const SQRT3: f32 = 1.732; // sqrt() is not a const fn
const TIME_STEP_SI: f32 = CELL_SIZE_SI / (SQRT3 * SOUND_SPEED_SI);
const RELAXATION_FACTOR: f32 =
    1.0 / (0.5 + (SQRT3 * KIN_VISCOSITY_SI) / (SOUND_SPEED_SI * CELL_SIZE_SI));
const CELLS: usize = PHYS_COLS * PHYS_ROWS;
const DEFAULT_VELOCITY: Vector2<f32> = Vector2::new(
    DEFAULT_VELOCITY_SI.x / (SQRT3 * SOUND_SPEED_SI),
    DEFAULT_VELOCITY_SI.y / (SQRT3 * SOUND_SPEED_SI),
);

// Lattice
#[rustfmt::skip]
const C_VELS: [Vector2<f32>; 9] = {
    const fn v(dr: isize, dc: isize) -> Vector2<f32> {
        Vector2::new(dr as f32, dc as f32)
    }
    [v(-1, -1), v(-1, 0), v(-1, 1),
     v( 0, -1), v( 0, 0), v( 0, 1),
     v( 1, -1), v( 1, 0), v( 1, 1)]
};
const WEIGHT: [f32; 9] = {
    let a = 1.0 / 9.0;
    let b = 1.0 / 36.0;
    [b, a, b, a, 4.0 / 9.0, a, b, a, b]
};

pub struct Physics {
    lattice: Box<[[[f32; 9]; PHYS_COLS]; PHYS_ROWS]>,
    blocked: Box<[[bool; PHYS_COLS]; PHYS_ROWS]>,
    // Measured in abstract time units. A cell has width 1 space unit and the
    // speed of sound is roughly 1 space unit per time unit.
    time_elapsed: usize,
}

impl Physics {
    pub fn is_blocked_at(&self, (r, c): (usize, usize)) -> bool {
        self.blocked[r][c]
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
            .fold(Vector2::zero(), |acc, (i, &val)| acc + C_VELS[i] * val)
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
        cells_iter_mut(&mut lattice).for_each(|cell| *cell = equiv);
        let blocked = vec![false; CELLS];
        let blocked: Box<[bool; CELLS]> = blocked.into_boxed_slice().try_into().unwrap();
        unsafe {
            Self {
                lattice,
                blocked: std::mem::transmute(blocked),
                time_elapsed: 0,
            }
        }
    }
    fn new_lattice() -> Box<[[[f32; 9]; PHYS_COLS]; PHYS_ROWS]> {
        let lattice = vec![0.0; 9 * CELLS];
        let lattice: Box<[f32; 9 * CELLS]> = lattice.into_boxed_slice().try_into().unwrap();
        unsafe { std::mem::transmute(lattice) }
    }
    pub fn lbm_collide(&mut self) {
        cells_iter_mut(&mut self.lattice).for_each(|cell| collide_cell(cell));
    }
    pub fn lbm_stream(&mut self) {
        let mut new_lattice = Self::new_lattice();
        for (r0, c0) in cells().filter(|&cell| !self.is_blocked_at(cell)) {
            for &(dr0, dc0) in directions() {
                let (r, c, dr, dc) = self.go((r0, c0), (dr0, dc0));
                assert!(r < PHYS_ROWS);
                assert!(c < PHYS_COLS);
                assert!(!self.is_blocked_at((r, c)));
                let old_vel = (3 * dr0 + dc0 + 4) as usize;
                let new_vel = (3 * dr + dc + 4) as usize;
                new_lattice[r][c][new_vel] += self.lattice[r0][c0][old_vel];
            }
        }
        for r in 0..PHYS_ROWS {
            set_vel(&mut new_lattice[r][0], DEFAULT_VELOCITY);
            set_vel(&mut new_lattice[r][PHYS_COLS - 1], DEFAULT_VELOCITY);
        }
        self.lattice = new_lattice;
        self.time_elapsed += 1;
    }
    fn go(&self, (r, c): (usize, usize), (dr, dc): (isize, isize)) -> (usize, usize, isize, isize) {
        // Bounce at upper and lower wall
        let (nr, dr) = {
            let nr = r as isize + dr;
            if nr == -1 {
                (0, 1)
            } else if nr == PHYS_ROWS as isize {
                (PHYS_ROWS - 1, -1)
            } else {
                (nr as usize, dr)
            }
        };
        // Periodic at left and right wall
        let nc = {
            let nc = c as isize + dc;
            if nc == -1 {
                PHYS_COLS - 1
            } else if nc == PHYS_COLS as isize {
                0
            } else {
                nc as usize
            }
        };
        // Bounce against blocked cells
        let blocked_both = self.is_blocked_at((nr, nc));
        let blocked_r = self.is_blocked_at((nr, c));
        let blocked_c = self.is_blocked_at((r, nc));
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
    pub fn get_velocity_vorticity_speed(
        &self,
        velocity: &mut [[Vector2<f32>; PHYS_COLS]; PHYS_ROWS],
        vorticity: &mut [[f32; PHYS_COLS]; PHYS_ROWS],
        speed: &mut [[f32; PHYS_COLS]; PHYS_ROWS],
    ) -> () {
        for (r, c) in cells() {
            let vel = &mut velocity[r][c];
            *vel = self.momentum_at((r, c)) / self.mass_at((r, c));
            if !vel.is_finite() {
                *vel = Vector2::zero();
            }
            speed[r][c] = vel.magnitude();
        }
        for (r, c) in cells() {
            if self.is_blocked_at((r, c)) {
                continue;
            }
            let mut res = 0.0;
            if r != 0 && r != PHYS_ROWS - 1 {
                res -= velocity[r - 1][c].y;
                res += velocity[r + 1][c].y;
            }
            if c != 0 && c != PHYS_COLS - 1 {
                res += velocity[r][c - 1].x;
                res -= velocity[r][c + 1].x;
            }
            vorticity[r][c] = if res.is_finite() { res } else { 0.0 };
        }
    }
}

fn set_vel(cell: &mut [f32; 9], vel: Vector2<f32>) {
    let mass: f32 = cell.iter().sum();
    cell.iter_mut()
        .zip(equilibrium(vel))
        .for_each(|(cell, eq)| *cell = mass * eq);
}

fn cells_iter_mut(
    lattice: &mut [[[f32; 9]; PHYS_COLS]; PHYS_ROWS],
) -> impl Iterator<Item = &mut [f32; 9]> {
    lattice
        .iter_mut()
        .map(|lattice_row| lattice_row.iter_mut())
        .flatten()
}

#[rustfmt::skip]
fn directions() -> &'static [(isize, isize)] {
    &[(-1, -1), (-1, 0), (-1, 1),
      ( 0, -1), ( 0, 0), ( 0, 1),
      ( 1, -1), ( 1, 0), ( 1, 1)]
}

pub fn cells() -> impl Iterator<Item = (usize, usize)> {
    (0..PHYS_ROWS)
        .map(|r| (0..PHYS_COLS).map(move |c| (r, c)))
        .flatten()
}

fn equilibrium(vel: Vector2<f32>) -> impl Iterator<Item = f32> {
    (0..9)
        .map(move |v| {
            let along = C_VELS[v].dot(vel);

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
    assert!(mass.is_finite());
    if mass.is_zero() {
        return;
    }
    assert!(mass.is_sign_positive());

    let momentum: Vector2<f32> = cell.iter().enumerate().map(|(i, &p)| C_VELS[i] * p).sum();
    assert!(momentum.is_finite());

    let vel = momentum / mass;
    assert!(vel.is_finite());

    cell.iter_mut()
        .zip(equilibrium(vel).map(|eq| mass * eq))
        .for_each(|(f, f_eq)| *f += RELAXATION_FACTOR * (f_eq - *f));
}
