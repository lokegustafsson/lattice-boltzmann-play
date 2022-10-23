use cgmath::Vector2;
use std::time::Duration;

pub struct Config {
    /// pixels length per simulation cell length, for rendering
    pub subcells: usize,
    /// physics updates for every render update
    pub substeps: usize,
    pub min_duration_between_logs: Duration,
    pub fps_limit: usize,
    pub physics: PhysicsConfig,
}
pub struct PhysicsConfig {
    pub grid_cols: usize,
    pub grid_rows: usize,
    pub si_cell_size: f32, // m
    pub mach_number_default_velocity: Vector2<f32>,
    pub reynolds_number_undisturbed_flow: f32,
}
pub const CONFIG: Config = Config {
    subcells: 4,
    substeps: 20,
    min_duration_between_logs: Duration::from_secs(1),
    fps_limit: 60,
    physics: PhysicsConfig {
        grid_cols: 400,
        grid_rows: 250,
        si_cell_size: 0.001,
        mach_number_default_velocity: Vector2::new(0.0, 0.1),
        reynolds_number_undisturbed_flow: 1e5,
    },
};

// FIXME: IS THIS EVEN TRUE?
// LBM restricts us to relatively low reynolds numbers. Also, since this implementation is very
// slow, we need to keep a significant fraction of the speed of sound to see any action. Hence,
// this simulation simulates a few grams of honey flowing through a tube the size of a water
// bottle at 100 m/s, in extreme slow motion
