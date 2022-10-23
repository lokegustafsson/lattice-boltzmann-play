use cgmath::Vector2;
use std::time::Duration;

pub struct Config {
    /// pixels length per simulation cell length, for rendering
    pub subcells: usize,
    /// physics updates for every render update
    pub substeps: usize,
    pub min_duration_between_logs: Duration,
    pub fps_limit: usize,
    pub test_particle_speed: f32,
    pub physics: PhysicsConfig,
}
pub struct PhysicsConfig {
    pub sound_speed_si: f32,   // m/s
    pub kin_viscosity_si: f32, // m^2/s
    pub density_si: f32,       // kg/m^3
    pub grid_cols: usize,
    pub grid_rows: usize,
    pub cell_size_si: f32, // m
    pub default_velocity_si: Vector2<f32>,
}
pub const CONFIG: Config = Config {
    subcells: 4,
    substeps: 20,
    min_duration_between_logs: Duration::from_secs(1),
    fps_limit: 600,
    test_particle_speed: 2.0,
    physics: PhysicsConfig {
        // Physical (honey, roughly?)
        sound_speed_si: 1400.0,
        kin_viscosity_si: 0.005,
        density_si: 1000.0,
        grid_cols: 400,
        grid_rows: 250,
        cell_size_si: 0.001,
        default_velocity_si: Vector2::new(0.0, 100.0),
    },
};

// FIXME: IS THIS EVEN TRUE?
// LBM restricts us to relatively low reynolds numbers. Also, since this implementation is very
// slow, we need to keep a significant fraction of the speed of sound to see any action. Hence,
// this simulation simulates a few grams of honey flowing through a tube the size of a water
// bottle at 100 m/s, in extreme slow motion
