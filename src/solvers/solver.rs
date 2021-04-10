///The Interface for wave propagation solvers

enum WetDryState {
    DryDry,
    WetWet,
    WetDryInundation,
    WetDryWall,
    WetDryWallInundation,
    DryWetInundation,
    DryWetWall,
    DryWetWallInundation,
}

pub struct SolverData {
    pub dry_tolerance: f32,
    pub g: f32,
    pub zero_tolerance: f32,

    pub h_left: f32,
    pub h_right: f32,
    pub hu_left: f32,
    pub hu_right: f32,
    pub b_left: f32,
    pub u_right: f32,
    pub u_left: f32,

    pub wet_dry_state: WetDryState,
}

pub struct Updates {
    h_update_left: f32,
    h_update_right: f32,
    hu_update_left: f32,
    hu_update_right: f32,
    max_wavespeed: f32,
}

impl Updates {
    pub fn new() -> Updates {
        Updates {
            h_update_left: 0.0,
            h_update_right: 0.0,
            hu_update_left: 0.0,
            hu_update_right: 0.0,
            max_wavespeed: 0.0,
        }
    }
}

trait Solver {
    fn compute_net_updates(&self) -> Updates;
}

impl SolverData {
    pub fn new(i_dry_tolerance: f32, i_gravity: f32, i_zero_tolerance: f32) -> SolverData {
        SolverData {
            dry_tolerance: i_dry_tolerance,
            g: i_gravity,
            zero_tolerance: i_zero_tolerance,
            h_left: 0.0,
            h_right: 0.0,
            hu_left: 0.0,
            hu_right: 0.0,
            b_left: 0.0,
            u_right: 0.0,
            u_left: 0.0,
        }
    }
}