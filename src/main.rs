mod framebuffer;
mod physics;
mod visual;

use crate::{
    framebuffer::FrameBuffer,
    physics::{Physics, PHYS_COLS, PHYS_ROWS},
    visual::TestParticle,
};
use cgmath::{Vector2, Zero};
use minifb::{Key, KeyRepeat, Window, WindowOptions};
use std::{
    ops::{AddAssign, Div, Mul},
    time::{Duration, Instant},
};

const PIXELS_PER_CELL: usize = 8;
const SUBSTEPS: usize = 10;
const MIN_DURATION_BETWEEN_LOGS: Duration = Duration::from_secs(1);
const FPS_LIMIT: usize = 60;

pub struct Sim {
    test_particles: Vec<TestParticle>,
    physics: Physics,
    ui: UiState,
    measured: Measured,
}
struct UiState {
    tick_index: usize,
    last_tick_instant: Instant,
    last_logged_instant: Instant,
    paint_cells_with_filled: bool,
    last_tick_had_mouse_down: bool,
    show_vorticity: bool,
    visualization_sensitivity: i8,
}
struct InputSnapshot {
    mouse_is_down: bool,
    mouse_pos: (usize, usize),
    pressed_space: bool,
    pressed_up: bool,
    pressed_down: bool,
}
struct Measured {
    velocity: Box<[[Vector2<f32>; PHYS_COLS]; PHYS_ROWS]>,
    vorticity: Box<[[f32; PHYS_COLS]; PHYS_ROWS]>,
    speed: Box<[[f32; PHYS_COLS]; PHYS_ROWS]>,
}
impl Sim {
    fn new() -> Self {
        Self {
            test_particles: (0..1000)
                .map(|_| TestParticle::new_random_anywhere())
                .collect(),
            physics: Physics::new(),
            ui: UiState {
                tick_index: 0,
                last_tick_instant: Instant::now(),
                last_logged_instant: Instant::now(),
                paint_cells_with_filled: true,
                last_tick_had_mouse_down: false,
                show_vorticity: false,
                visualization_sensitivity: 1,
            },
            measured: Measured {
                velocity: Box::new([[Vector2::zero(); PHYS_COLS]; PHYS_ROWS]),
                vorticity: Box::new([[0.0; PHYS_COLS]; PHYS_ROWS]),
                speed: Box::new([[0.0; PHYS_COLS]; PHYS_ROWS]),
            },
        }
    }
    fn handle_input(
        &mut self,
        InputSnapshot {
            mouse_is_down,
            mouse_pos: (wx, wy),
            pressed_space,
            pressed_up,
            pressed_down,
        }: InputSnapshot,
    ) {
        if mouse_is_down {
            self.physics.update_cell_blocked_status(
                FloatingCoords::from_window((wy, wx)).to_physics_cell(),
                self.ui.paint_cells_with_filled,
            );
        }
        if self.ui.last_tick_had_mouse_down && !mouse_is_down {
            self.ui.paint_cells_with_filled ^= true;
        }
        self.ui.show_vorticity ^= pressed_space;
        self.ui.visualization_sensitivity += pressed_up as i8;
        self.ui.visualization_sensitivity -= pressed_down as i8;
    }
    fn maybe_print_log(&mut self) {
        self.ui.tick_index += 1;
        let now = Instant::now();
        if now.duration_since(self.ui.last_logged_instant) > MIN_DURATION_BETWEEN_LOGS {
            println!(
                "Mass: {:.2} gram, Elapsed: {:.2} milliseconds, {:.1} FPS, MLUPS: {:.2}",
                self.physics.total_mass() * 1e3,
                self.physics.time_elapsed() * 1e3,
                1.0 / Instant::now()
                    .duration_since(self.ui.last_tick_instant)
                    .as_secs_f32(),
                (PHYS_ROWS * PHYS_COLS * SUBSTEPS) as f32
                    / (1e6
                        * Instant::now()
                            .duration_since(self.ui.last_tick_instant)
                            .as_secs_f32()),
            );
            self.ui.last_logged_instant = now;
        }
        self.ui.last_tick_instant = now;
    }
}

fn main() {
    let mut framebuffer = FrameBuffer::new();
    let mut window = Window::new(
        "A basic Lattice Boltzmann simulation",
        FrameBuffer::WIDTH,
        FrameBuffer::HEIGHT,
        WindowOptions::default(),
    )
    .unwrap();
    window.limit_update_rate(Some(std::time::Duration::from_nanos(
        1_000_000_000 / FPS_LIMIT as u64 / SUBSTEPS as u64,
    )));

    let mut sim = Sim::new();
    println!("Initial mass: {:.2} gram", sim.physics.total_mass() * 1e3);

    while window.is_open() && !window.is_key_down(Key::Escape) {
        for i in 0..SUBSTEPS {
            sim.handle_input(InputSnapshot {
                mouse_is_down: window.get_mouse_down(minifb::MouseButton::Left),
                mouse_pos: {
                    let (fx, fy) = window.get_mouse_pos(minifb::MouseMode::Clamp).unwrap();
                    (fx as usize, fy as usize)
                },
                pressed_space: window.is_key_pressed(Key::Space, KeyRepeat::No),
                pressed_up: window.is_key_pressed(Key::Up, KeyRepeat::No),
                pressed_down: window.is_key_pressed(Key::Down, KeyRepeat::No),
            });
            if i != SUBSTEPS - 1 {
                window.update();
            }
            sim.physics.lbm_collide();
            sim.physics.lbm_stream();
        }
        sim.maybe_print_log();
        sim.physics.get_velocity_vorticity_speed(
            &mut sim.measured.velocity,
            &mut sim.measured.vorticity,
            &mut sim.measured.speed,
        );

        framebuffer.update(
            sim.physics.blocked(),
            if sim.ui.show_vorticity {
                &sim.measured.vorticity
            } else {
                &sim.measured.speed
            },
            sim.ui.visualization_sensitivity as f32,
        );
        {
            let vel = &sim.measured.velocity;
            sim.test_particles.iter_mut().for_each(|t| t.advance(vel));
            for test_particle in &sim.test_particles {
                test_particle.paint_on(&mut framebuffer, sim.ui.show_vorticity);
            }
            for test_particle in &mut sim.test_particles {
                if test_particle.is_outside_bounds() {
                    *test_particle = TestParticle::new_random_inflow();
                }
            }
        }

        window
            .update_with_buffer(
                framebuffer.as_flat(),
                FrameBuffer::WIDTH,
                FrameBuffer::HEIGHT,
            )
            .unwrap();
    }
}

#[derive(Clone, Copy)]
pub struct FloatingCoords {
    r: f32,
    c: f32,
}
impl FloatingCoords {
    fn from_window((r, c): (usize, usize)) -> Self {
        Self {
            r: r as f32 / (FrameBuffer::HEIGHT - 1) as f32,
            c: c as f32 / (FrameBuffer::WIDTH - 1) as f32,
        }
    }
    fn to_window(self) -> (usize, usize) {
        let r = (self.r * (FrameBuffer::HEIGHT - 1) as f32) as usize;
        let c = (self.c * (FrameBuffer::WIDTH - 1) as f32) as usize;
        (r, c)
    }
    fn to_physics_cell(self) -> (usize, usize) {
        let ri = ((self.r * (PHYS_ROWS - 1) as f32) + 0.5) as usize;
        let ci = ((self.c * (PHYS_COLS - 1) as f32) + 0.5) as usize;
        (ri, ci)
    }

    fn sample_field_interpolated<T>(self, field: &[[T; PHYS_COLS]; PHYS_ROWS]) -> T
    where
        T: AddAssign<T> + Mul<f32, Output = T> + Div<f32, Output = T>,
        T: Zero + Copy,
    {
        let r = self.r * (PHYS_ROWS - 1) as f32;
        let c = self.c * (PHYS_COLS - 1) as f32;
        let ri = (r + 0.5) as usize;
        let ci = (c + 0.5) as usize;

        let mut res: T = T::zero();
        let mut tot_weight: f32 = 0.0;
        for dr in 0..3 {
            for dc in 0..3 {
                let r_point = (ri + dr).wrapping_sub(1);
                let c_point = (ci + dc).wrapping_sub(1);
                if r_point < PHYS_ROWS && c_point < PHYS_COLS {
                    let r_delta = r_point as f32 - r;
                    let c_delta = c_point as f32 - c;
                    let weight: f32 = (1.0 - r_delta * r_delta - c_delta * c_delta).max(0.0);

                    res += field[r_point][c_point] * weight;
                    tot_weight += weight;
                }
            }
        }
        res / tot_weight
    }
}
