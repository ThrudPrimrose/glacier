use std::f32;

///The Interface for wave propagation solvers
pub mod fwave;

/*
pub enum WetDryState {
    DryDry,
    WetWet,
    WetDryInundation,
    WetDryWall,
    WetDryWallInundation,
    DryWetInundation,
    DryWetWall,
    DryWetWallInundation,
}
*/

pub struct SolverData {
    pub dry_tolerance: f32,
    pub g: f32,
    pub half_g: f32,
    pub sqrt_g: f32,
    pub zero_tolerance: f32,
    /*
    pub h_left: f32,
    pub h_right: f32,
    pub hu_left: f32,
    pub hu_right: f32,
    pub b_left: f32,
    pub b_right: f32,
    pub u_right: f32,
    pub u_left: f32,

    pub wet_dry_state: WetDryState,
    */
}

#[derive(Clone, Copy)]
pub struct Update {
    pub h_update_left: f32,
    pub h_update_right: f32,
    pub hu_update_left: f32,
    pub hu_update_right: f32,
    pub max_wavespeed: f32,
}

impl Update {
    pub fn new() -> Update {
        Update {
            h_update_left: 0.0,
            h_update_right: 0.0,
            hu_update_left: 0.0,
            hu_update_right: 0.0,
            max_wavespeed: 0.0,
        }
    }

    pub fn wavespeed_only(wavespeed: f32) -> Update {
        Update {
            h_update_left: 0.0,
            h_update_right: 0.0,
            hu_update_left: 0.0,
            hu_update_right: 0.0,
            max_wavespeed: wavespeed,
        }
    }
}

impl SolverData {
    pub fn new(i_dry_tolerance: f32, i_gravity: f32, i_zero_tolerance: f32) -> SolverData {
        SolverData {
            dry_tolerance: i_dry_tolerance,
            g: i_gravity,
            half_g: i_gravity / 2.0,
            sqrt_g: i_gravity.sqrt(),
            zero_tolerance: i_zero_tolerance,
            /*
            h_left: 0.0,
            h_right: 0.0,
            hu_left: 0.0,
            hu_right: 0.0,
            b_left: 0.0,
            b_right: 0.0,
            u_right: 0.0,
            u_left: 0.0,
            wet_dry_state: WetDryState::WetWet,
            */
        }
    }

    pub fn default() -> SolverData {
        let g: f32 = 9.81;
        SolverData {
            dry_tolerance: 0.1,
            g: g,
            half_g: g / 2.0,
            sqrt_g: g.sqrt(),
            zero_tolerance: 0.0000001,
            /*
            h_left: 0.0,
            h_right: 0.0,
            hu_left: 0.0,
            hu_right: 0.0,
            b_left: 0.0,
            b_right: 0.0,
            u_right: 0.0,
            u_left: 0.0,
            wet_dry_state: WetDryState::WetWet,
            */
        }
    }
}

pub trait Solver {
    /*
    fn set_values(
        &mut self,
        h_left: f32,
        h_right: f32,
        hu_left: f32,
        hu_right: f32,
        b_left: f32,
        b_right: f32,
    ) -> ();
    */
    fn compute_net_updates(
        &self,
        h_left: f32,
        h_right: f32,
        hu_right: f32,
        hu_left: f32,
        b_left: f32,
        b_right: f32,
    ) -> Update;
}
