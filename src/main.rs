mod lattice;
mod visual;

use lattice::{Lattice, COLS, ROWS};
use minifb::{Key, KeyRepeat, Window, WindowOptions};
use rayon::prelude::*;
use visual::{float_to_sim, interpolate, speed, two_color, velocity, vorticity, Trace};

const WIDTH: usize = 10 * COLS;
const HEIGHT: usize = 10 * ROWS;

fn main() {
    let mut buffer: Vec<u32> = vec![0; WIDTH * HEIGHT];

    let mut window = Window::new(
        "A basic Lattice Boltzmann simulation",
        WIDTH,
        HEIGHT,
        WindowOptions::default(),
    )
    .unwrap_or_else(|e| {
        panic!("{}", e);
    });

    // Limit to max ~60 fps update rate
    window.limit_update_rate(Some(std::time::Duration::from_micros(16600)));

    let mut traces: Vec<Trace> = (0..1000).map(|_| Trace::initial()).collect();
    let mut lattice = Lattice::new();
    lattice
        .iter_mut()
        .for_each(|cell| cell[4] = 1.0 + 0.1 * rand::random::<f32>());
    println!("Initial mass: {}", lattice.total_mass());

    let mut last = String::new();
    let mut unblocking = true;
    let mut clicking = false;
    let mut show_vorticity = false;

    while window.is_open() && !window.is_key_down(Key::Escape) {
        if !clicking && window.get_mouse_down(minifb::MouseButton::Left) {
            clicking = true;
            let (wx, wy) = window.get_mouse_pos(minifb::MouseMode::Clamp).unwrap();
            let cell = float_to_sim(window_to_float((wy as usize, wx as usize)));
            unblocking = lattice.blocked(cell);
        } else if clicking && !window.get_mouse_down(minifb::MouseButton::Left) {
            clicking = false;
        }
        if clicking {
            let (wx, wy) = window.get_mouse_pos(minifb::MouseMode::Clamp).unwrap();
            let cell = float_to_sim(window_to_float((wy as usize, wx as usize)));
            if lattice.blocked(cell) == unblocking {
                lattice.toggle_block(cell);
            }
        }
        if window.is_key_pressed(Key::Space, KeyRepeat::No) {
            show_vorticity ^= true;
        }

        lattice.collide();
        lattice.stream();
        let mass = lattice.total_mass();
        let msg = format!("{:.2e}", mass);
        if msg != last {
            last = msg;
            println!("Mass: {}", last);
        }

        let vorticity = vorticity(&lattice);
        let velocity = velocity(&lattice);
        let speed = speed(&velocity);

        pixels()
            .into_par_iter()
            .map(|win_pos| {
                let float_pos = window_to_float(win_pos);
                let cell = float_to_sim(float_pos);
                if lattice.blocked(cell) {
                    return (win_pos, 0x00FF00);
                }

                let field = if show_vorticity { &vorticity } else { &speed };
                let value = interpolate(field, float_pos);
                const CUTOFF: f32 = 0.00001;
                let res = value.signum()
                    * (((value.abs() + CUTOFF).ln() - CUTOFF.ln()) / (-CUTOFF.ln())).min(1.0);
                return (win_pos, two_color(res));
            })
            .collect::<Vec<_>>()
            .iter()
            .for_each(|((r, c), color)| buffer[WIDTH * r + c] = *color);

        traces.iter_mut().for_each(|t| t.advance(&velocity));
        for trace in &traces {
            let win_pos = float_to_window(trace.get());
            for dr in 0..3 {
                for dc in 0..3 {
                    let r = win_pos.0 + dr - 1;
                    let c = win_pos.1 + dc - 1;
                    if r < HEIGHT && c < WIDTH {
                        buffer[WIDTH * r + c] = if show_vorticity { 0x000000 } else { 0xFF0000 };
                    }
                }
            }
        }
        for trace in &mut traces {
            if !trace.inside() {
                *trace = Trace::new();
            }
        }

        window.update_with_buffer(&buffer, WIDTH, HEIGHT).unwrap();
    }
}

fn pixels() -> impl ParallelIterator<Item = (usize, usize)> {
    (0..HEIGHT)
        .into_par_iter()
        .map(|r| (0..WIDTH).into_par_iter().map(move |c| (r, c)))
        .flatten()
}

fn window_to_float((r, c): (usize, usize)) -> (f32, f32) {
    let rf = r as f32 / (HEIGHT - 1) as f32;
    let cf = c as f32 / (WIDTH - 1) as f32;
    (rf, cf)
}

fn float_to_window((rf, cf): (f32, f32)) -> (usize, usize) {
    let r = (rf * (HEIGHT - 1) as f32) as usize;
    let c = (cf * (WIDTH - 1) as f32) as usize;
    (r, c)
}
