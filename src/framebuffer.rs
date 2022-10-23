use crate::config::CONFIG as C;
use rayon::prelude::*;

const GRID_COLS: usize = C.physics.grid_cols;
const GRID_ROWS: usize = C.physics.grid_rows;
type NestedArray = [[[[u32; C.subcells]; GRID_COLS]; C.subcells]; GRID_ROWS];
pub struct FrameBuffer {
    inner: Box<NestedArray>,
}
impl FrameBuffer {
    pub const WIDTH: usize = C.subcells * GRID_COLS;
    pub const HEIGHT: usize = C.subcells * GRID_ROWS;

    const SIZE: usize = C.subcells * GRID_COLS * C.subcells * GRID_ROWS;

    pub fn new() -> Self {
        Self {
            inner: Box::new([[[[0; C.subcells]; GRID_COLS]; C.subcells]; GRID_ROWS]),
        }
    }
    pub fn as_flat(&self) -> &[u32; Self::SIZE] {
        unsafe { &*(self.inner.as_ref() as *const NestedArray as *const [u32; Self::SIZE]) }
    }
    pub fn set_pixel(&mut self, r: usize, c: usize, color: u32) {
        if r > Self::HEIGHT || c > Self::WIDTH {
            dbg!(r, c, Self::HEIGHT, Self::WIDTH);
        }
        self.inner[r / C.subcells][r % C.subcells][c / C.subcells][c % C.subcells] = color
    }
    pub fn update(
        &mut self,
        blocked: &[[bool; GRID_COLS]; GRID_ROWS],
        field: &[[f32; GRID_COLS]; GRID_ROWS],
        visualization_sensitivity: f32,
    ) {
        self.inner
            .par_iter_mut()
            .enumerate()
            .for_each(|(r, block_row)| {
                for (cell_r, row) in block_row.iter_mut().enumerate() {
                    for (c, block_col) in row.iter_mut().enumerate() {
                        for (cell_c, cell) in block_col.iter_mut().enumerate() {
                            let color = if blocked[r][c] {
                                0x00FF00
                            } else {
                                let val = sample_field(field, r, c, cell_r, cell_c)
                                    * visualization_sensitivity
                                    * 255.0;
                                let intensity = (val.abs() as u32).min(255);
                                0xFFFFFF - (if val < 0.0 { 0x000101 } else { 0x010100 }) * intensity
                            };
                            *cell = color;
                        }
                    }
                }
            });
    }
}

fn sample_field(
    field: &[[f32; GRID_COLS]; GRID_ROWS],
    r: usize,
    c: usize,
    cell_r: usize,
    cell_c: usize,
) -> f32 {
    const HALFCELL: usize = C.subcells / 2;
    const TC: f32 = C.subcells as f32;
    if r == 0 || r + 1 == GRID_ROWS || c == 0 || c + 1 == GRID_COLS {
        field[r][c]
    } else {
        match (cell_r >= HALFCELL, cell_c >= HALFCELL) {
            (false, false) => {
                let r2 = 0.5 + 0.5 / TC + (cell_r as f32) / TC;
                let c2 = 0.5 + 0.5 / TC + (cell_c as f32) / TC;
                let r1 = 1.0 - r2;
                let c1 = 1.0 - c2;
                r1 * (c1 * field[r - 1][c - 1] + c2 * field[r - 1][c])
                    + r2 * (c1 * field[r][c - 1] + c2 * field[r][c])
            }
            (false, true) => {
                let r2 = 0.5 + 0.5 / TC + (cell_r as f32) / TC;
                let c2 = -0.5 + 0.5 / TC + (cell_c as f32) / TC;
                let r1 = 1.0 - r2;
                let c1 = 1.0 - c2;
                r1 * (c1 * field[r - 1][c] + c2 * field[r - 1][c + 1])
                    + r2 * (c1 * field[r][c] + c2 * field[r][c + 1])
            }
            (true, false) => {
                let r2 = -0.5 + 0.5 / TC + (cell_r as f32) / TC;
                let c2 = 0.5 + 0.5 / TC + (cell_c as f32) / TC;
                let r1 = 1.0 - r2;
                let c1 = 1.0 - c2;
                r1 * (c1 * field[r][c - 1] + c2 * field[r][c])
                    + r2 * (c1 * field[r + 1][c - 1] + c2 * field[r + 1][c])
            }
            (true, true) => {
                let r2 = -0.5 + 0.5 / TC + (cell_r as f32) / TC;
                let c2 = -0.5 + 0.5 / TC + (cell_c as f32) / TC;
                let r1 = 1.0 - r2;
                let c1 = 1.0 - c2;
                r1 * (c1 * field[r][c] + c2 * field[r][c + 1])
                    + r2 * (c1 * field[r + 1][c] + c2 * field[r + 1][c + 1])
            }
        }
    }
}
