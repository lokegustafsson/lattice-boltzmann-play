use crate::lattice::{cells, Lattice, COLS, ROWS};
use cgmath::{prelude::*, Vector2, Zero};
use std::convert::TryInto;
use std::ops::{AddAssign, Div, Mul};

pub fn float_to_sim((r, c): (f32, f32)) -> (usize, usize) {
    let ri = ((r * (ROWS - 1) as f32) + 0.5) as usize;
    let ci = ((c * (COLS - 1) as f32) + 0.5) as usize;
    (ri, ci)
}

pub fn two_color(value: f32) -> u32 {
    let intensity = ((1.0 - 2.0 * (value.abs() - 0.5).abs()) * 255.0) as u32;
    if value < -0.5 {
        0x010000 * intensity
    } else if value < 0.0 {
        0xFFFFFF - 0x000101 * intensity
    } else if value < 0.5 {
        0xFFFFFF - 0x010100 * intensity
    } else {
        0x000001 * intensity
    }
}

pub fn interpolate<T>(field: &[[T; COLS]; ROWS], (r, c): (f32, f32)) -> T
where
    T: AddAssign<T> + Mul<f32, Output = T> + Div<f32, Output = T>,
    T: Zero + Copy,
{
    let r = r * (ROWS - 1) as f32;
    let c = c * (COLS - 1) as f32;
    let ri = (r + 0.5) as usize;
    let ci = (c + 0.5) as usize;

    let mut res: T = T::zero();
    let mut tot_weight: f32 = 0.0;
    for dr in 0..3 {
        for dc in 0..3 {
            let r_point = ri + dr - 1;
            let c_point = ci + dc - 1;
            if r_point < ROWS && c_point < COLS {
                let r_delta = r_point as f32 - r;
                let c_delta = c_point as f32 - c;
                let weight: f32 = (1.0 - r_delta.hypot(c_delta)).max(0.0);

                res += field[r_point][c_point] * weight;
                tot_weight += weight;
            }
        }
    }
    res / tot_weight
}

pub fn vorticity(lattice: &Lattice) -> Box<[[f32; COLS]; ROWS]> {
    let mut field: Box<[[f32; COLS]; ROWS]> = {
        let field = vec![0.0; ROWS * COLS];
        let field: Box<[f32; ROWS * COLS]> = field.into_boxed_slice().try_into().unwrap();
        unsafe { std::mem::transmute(field) }
    };
    for (r, c) in cells() {
        if lattice.blocked((r, c)) {
            continue;
        }
        let mut res = 0.0;
        if r != 0 && r != ROWS - 1 {
            res -= lattice.momentum_at((r - 1, c)).y;
            res += lattice.momentum_at((r + 1, c)).y;
        }
        if c != 0 && c != COLS - 1 {
            res += lattice.momentum_at((r, c - 1)).x;
            res -= lattice.momentum_at((r, c + 1)).x;
        }
        field[r][c] = if res.is_finite() { res } else { 0.0 };
    }
    field
}

pub fn velocity(lattice: &Lattice) -> Box<[[Vector2<f32>; COLS]; ROWS]> {
    let mut field: Box<[[Vector2<f32>; COLS]; ROWS]> = {
        let field = vec![Vector2::zero(); ROWS * COLS];
        let field: Box<[Vector2<f32>; ROWS * COLS]> = field.into_boxed_slice().try_into().unwrap();
        unsafe { std::mem::transmute(field) }
    };
    for cell in cells() {
        field[cell.0][cell.1] = lattice.momentum_at(cell) / lattice.mass_at(cell);
        if !field[cell.0][cell.1].is_finite() {
            field[cell.0][cell.1] = Vector2::zero();
        }
    }
    field
}

pub fn speed(velocity: &[[Vector2<f32>; COLS]; ROWS]) -> Box<[[f32; COLS]; ROWS]> {
    let field: Vec<f32> = velocity
        .iter()
        .map(|s| s.iter())
        .flatten()
        .map(|v| v.magnitude())
        .collect();
    let field: Box<[f32; ROWS * COLS]> = field.into_boxed_slice().try_into().unwrap();
    unsafe { std::mem::transmute(field) }
}

pub struct Trace(f32, f32);

impl Trace {
    pub fn initial() -> Self {
        Self(rand::random::<f32>(), rand::random::<f32>())
    }
    pub fn new() -> Self {
        Self(rand::random::<f32>(), 0.01)
    }
    pub fn get(&self) -> (f32, f32) {
        (self.0, self.1)
    }
    pub fn inside(&self) -> bool {
        0.0 <= self.0 && self.0 <= 1.0 && 0.01 <= self.1 && self.1 <= 0.99
    }
    pub fn advance(&mut self, velocity: &[[Vector2<f32>; COLS]; ROWS]) {
        let vel = interpolate(velocity, (self.0, self.1));
        const SPEED: f32 = 3.0;
        self.0 += SPEED * vel.x / ROWS as f32;
        self.1 += SPEED * vel.y / COLS as f32;
    }
}
