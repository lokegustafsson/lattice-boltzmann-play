use crate::{
    framebuffer::FrameBuffer,
    physics::{PHYS_COLS, PHYS_ROWS},
    FloatingCoords, Sim,
};
use cgmath::Vector2;

pub struct TestParticle {
    pos: FloatingCoords,
}

impl TestParticle {
    pub fn new_random_anywhere() -> Self {
        Self {
            pos: FloatingCoords {
                r: rand::random::<f32>(),
                c: rand::random::<f32>(),
            },
        }
    }
    pub fn new_random_inflow() -> Self {
        Self {
            pos: FloatingCoords {
                r: rand::random::<f32>(),
                c: 0.01,
            },
        }
    }
    pub fn get(&self) -> FloatingCoords {
        self.pos
    }
    pub fn is_outside_bounds(&self) -> bool {
        !(0.0 <= self.pos.r && self.pos.r <= 1.0 && 0.01 <= self.pos.c && self.pos.c <= 0.99)
    }
    pub fn paint_on(&self, framebuffer: &mut FrameBuffer, show_vorticity: bool) {
        let (wr, wc) = self.get().to_window();
        for dr in 0..3 {
            for dc in 0..3 {
                let r = (wr + dr).wrapping_sub(1);
                let c = (wc + dc).wrapping_sub(1);
                if r < FrameBuffer::HEIGHT && c < FrameBuffer::WIDTH {
                    framebuffer.set_pixel(r, c, if show_vorticity { 0x000000 } else { 0xFF0000 });
                }
            }
        }
    }
    pub fn advance(&mut self, velocity: &[[Vector2<f32>; PHYS_COLS]; PHYS_ROWS]) {
        let vel = self.pos.sample_field_interpolated(velocity);
        const SPEED: f32 = 2.0;
        self.pos.r += SPEED * vel.x / PHYS_ROWS as f32;
        self.pos.c += SPEED * vel.y / PHYS_COLS as f32;
    }
}
