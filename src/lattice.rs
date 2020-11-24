use cgmath::{prelude::*, Vector2};
use std::convert::TryInto;
use std::iter::Iterator;

pub const COLS: usize = 100;
pub const ROWS: usize = 60;
const CELLS: usize = COLS * ROWS;

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

pub struct Lattice {
    lattice: Box<[[[f32; 9]; COLS]; ROWS]>,
    blocked: Box<[[bool; COLS]; ROWS]>,
}

impl Lattice {
    pub fn toggle_block(&mut self, (r, c): (usize, usize)) {
        self.blocked[r][c] ^= true;
    }
    pub fn blocked(&self, (r, c): (usize, usize)) -> bool {
        self.blocked[r][c]
    }
    pub fn total_mass(&self) -> f32 {
        let mut ans = 0.0;
        for (r, c) in cells().filter(|&c| !self.blocked(c)) {
            for v in 0..9 {
                ans += self.lattice[r][c][v];
            }
        }
        ans
    }
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut [f32; 9]> {
        self.lattice
            .iter_mut()
            .zip(self.blocked.iter())
            .map(|(lattice_row, blocked_row)| lattice_row.iter_mut().zip(blocked_row.iter()))
            .flatten()
            .filter(|(_, blocked)| !*blocked)
            .map(|(cell, _)| cell)
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

    pub fn new() -> Lattice {
        let blocked = vec![false; CELLS];
        let blocked: Box<[bool; CELLS]> = blocked.into_boxed_slice().try_into().unwrap();
        unsafe {
            Lattice {
                lattice: Lattice::new_lattice(),
                blocked: std::mem::transmute(blocked),
            }
        }
    }
    fn new_lattice() -> Box<[[[f32; 9]; COLS]; ROWS]> {
        let lattice = vec![0.0; 9 * CELLS];
        let lattice: Box<[f32; 9 * CELLS]> = lattice.into_boxed_slice().try_into().unwrap();
        unsafe { std::mem::transmute(lattice) }
    }
    pub fn collide(&mut self) {
        self.iter_mut().for_each(|cell| collide_cell(cell));
    }
    pub fn stream(&mut self) {
        let mut new_lattice = Lattice::new_lattice();
        for origin in cells().filter(|&cell| !self.blocked(cell)) {
            for jump in directions() {
                if let Some((nr, nc, dr, dc)) = self.go(origin, *jump) {
                    assert!(nr < ROWS);
                    assert!(nc < COLS);
                    assert!(!self.blocked((nr, nc)));
                    let vel = (3 * dr + dc + 4) as usize;
                    let amount =
                        self.lattice[origin.0][origin.1][(3 * jump.0 + jump.1 + 4) as usize];
                    new_lattice[nr][nc][vel] += amount;
                }
            }
        }
        for r in 0..ROWS {
            for v in 0..9 {
                new_lattice[r][0][v] = 0.0;
                new_lattice[r][COLS - 1][v] = 0.0;
            }
            new_lattice[r][0][4] = 10.0;
            new_lattice[r][COLS - 1][4] = 1.0;
        }
        self.lattice = new_lattice;
    }
    fn go(
        &self,
        (r, c): (usize, usize),
        (mut dr, dc): (isize, isize),
    ) -> Option<(usize, usize, isize, isize)> {
        let nc = c + dc as usize;
        if nc >= COLS {
            None
        } else {
            let nr = {
                let nr = r as isize + dr;
                if nr == -1 {
                    dr = 1;
                    0
                } else if nr == ROWS as isize {
                    dr = -1;
                    ROWS - 1
                } else {
                    nr as usize
                }
            };
            let blocked_both = self.blocked((nr, nc));
            let blocked_r = self.blocked((nr, c));
            let blocked_c = self.blocked((r, nc));
            if (blocked_r && blocked_c) || (!blocked_r && !blocked_c && blocked_both) {
                Some((r, c, -dr, -dc))
            } else if !blocked_both {
                Some((nr, nc, dr, dc))
            } else if blocked_c {
                Some((nr, c, dr, -dc))
            } else {
                Some((r, nc, -dr, dc))
            }
        }
    }
}

#[rustfmt::skip]
fn directions() -> &'static [(isize, isize)] {
    &[(-1, -1), (-1, 0), (-1, 1),
      ( 0, -1), ( 0, 0), ( 0, 1),
      ( 1, -1), ( 1, 0), ( 1, 1)]
}

pub fn cells() -> impl Iterator<Item = (usize, usize)> {
    (0..ROWS).map(|r| (0..COLS).map(move |c| (r, c))).flatten()
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

    let time_scale = 5.0;
    let sound2: f32 = (0..9).map(|i| WEIGHT[i] * C_VELS[i].magnitude2()).sum();

    for v in 0..9 {
        let uc = C_VELS[v].dot(vel);
        assert!(uc.is_finite());

        let uc_by_cs2 = uc / sound2;
        let equilibrium = WEIGHT[v]
            * mass
            * (1.0 + uc_by_cs2 + uc_by_cs2.powi(2) / 2.0 - vel.magnitude2() / sound2 / 2.0);
        assert!(equilibrium.is_finite());
        if equilibrium.is_sign_negative() {
            panic!("This breaks down when approaching/exceeding Mach 1... vel: {:?}, sound: {}", vel, sound2.sqrt());
        }

        cell[v] += (equilibrium - cell[v]) / time_scale;
    }
}
