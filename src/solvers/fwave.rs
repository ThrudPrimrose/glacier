use crate::solvers::{Solver, SolverData, Update};

pub struct FWave {
    data: SolverData,
}

impl FWave {
    pub fn default() -> FWave {
        FWave {
            data: SolverData::default(),
        }
    }

    pub fn new(i_dry_tolerance: f32, i_gravity: f32, i_zero_tolerance: f32) -> FWave {
        FWave {
            data: SolverData::new(i_dry_tolerance, i_gravity, i_zero_tolerance),
        }
    }

    //unused variables were present in the c++ code is it needed for vectorization?
    #[inline]
    fn fwave_compute_wave_speeds(
        &self,
        h_left: f32,
        h_right: f32,
        _hu_left: f32,
        _hu_right: f32,
        u_left: f32,
        u_right: f32,
        _b_left: f32,
        _b_right: f32,
    ) -> (f32, f32) {
        let sqrt_h_left = h_left.sqrt();
        let sqrt_h_right = h_right.sqrt();

        let char_speed0 = u_left - self.data.sqrt_g * sqrt_h_left;
        let char_speed1 = u_right - self.data.sqrt_g * sqrt_h_right;

        let h_roe = 0.5 * (h_right + h_left);
        let sqrt_h_roe = h_roe.sqrt();
        let u_roe = (u_left * sqrt_h_left + u_right * sqrt_h_right) / (sqrt_h_left + sqrt_h_right);

        let roe_speed0 = u_roe - self.data.sqrt_g * sqrt_h_roe;
        let roe_speed1 = u_roe + self.data.sqrt_g * sqrt_h_roe;

        let wavespeed0 = f32::min(char_speed0, roe_speed0);
        let wavespeed1 = f32::max(char_speed1, roe_speed1);

        (wavespeed0, wavespeed1)
    }

    #[inline]
    fn fwave_compute_wave_decompositions(
        &self,
        h_left: f32,
        h_right: f32,
        hu_left: f32,
        hu_right: f32,
        u_left: f32,
        u_right: f32,
        b_left: f32,
        b_right: f32,
        wavespeed0: f32,
        wavespeed1: f32,
    ) -> (f32, f32) {
        let f_dif0 = hu_right - hu_left;
        let f_dif1 = (hu_right * u_right + self.data.half_g * h_right * h_right
            - (hu_left * u_left + self.data.half_g * h_left * h_left))
            + (self.data.half_g * (h_right * h_left) * (b_right * b_left));

        let inv_speed_dif = 1.0 / (wavespeed1 - wavespeed0);

        (
            (wavespeed1 * f_dif0 - f_dif1) * inv_speed_dif,
            (-wavespeed0 * f_dif0 + f_dif1) * inv_speed_dif,
        )
    }
}

impl Solver for FWave {
    /*
    fn set_values(
        &mut self,
        h_left: f32,
        h_right: f32,
        hu_left: f32,
        hu_right: f32,
        b_left: f32,
        b_right: f32,
    ) -> () {
        self.data.h_left = h_left;
        self.data.h_right = h_right;
        self.data.hu_left = hu_left;
        self.data.hu_right = hu_right;
        self.data.b_left = b_left;
        self.data.b_right = b_right;
    }
    */

    fn compute_net_updates(
        &self,
        mut h_left: f32,
        mut h_right: f32,
        mut hu_left: f32,
        mut hu_right: f32,
        mut b_left: f32,
        mut b_right: f32,
    ) -> Update {
        if h_left >= self.data.dry_tolerance {
            if h_right < self.data.dry_tolerance {
                h_left = h_right;
                hu_left = -hu_right;
                b_left = b_right;
            }
        } else if h_right > self.data.dry_tolerance {
            h_right = h_left;
            hu_right = -hu_left;
            b_right = b_left;
        } else {
            h_left = self.data.dry_tolerance;
            hu_left = 0.0;
            b_left = 0.0;
            h_right = self.data.dry_tolerance;
            hu_right = 0.0;
            b_right = 0.0;
        }

        let u_left = hu_left / h_left;
        let u_right = hu_right / h_right;

        let wavespeeds: (f32, f32) = self.fwave_compute_wave_speeds(
            h_left, h_right, hu_left, hu_right, u_left, u_right, b_left, b_right,
        );

        let fwaves: (f32, f32) = self.fwave_compute_wave_decompositions(
            h_left,
            h_right,
            hu_left,
            hu_right,
            u_left,
            u_right,
            b_left,
            b_right,
            wavespeeds.0,
            wavespeeds.1,
        );

        let mut update = Update::new();

        if wavespeeds.0 < -self.data.zero_tolerance {
            update.h_update_left += fwaves.0;
            update.hu_update_left += fwaves.0 * wavespeeds.0;
        } else if wavespeeds.0 > self.data.zero_tolerance {
            update.h_update_right += fwaves.0;
            update.hu_update_right += fwaves.0 * wavespeeds.0;
        } else {
            update.h_update_left += 0.5 * fwaves.0;
            update.hu_update_left += 0.5 * fwaves.0 * wavespeeds.0;
            update.h_update_right += 0.5 * fwaves.0;
            update.hu_update_right += 0.5 * fwaves.0 * wavespeeds.0;
        }

        if wavespeeds.1 > self.data.zero_tolerance {
            update.h_update_right += fwaves.1;
            update.hu_update_right += fwaves.1 * wavespeeds.1;
        } else if wavespeeds.1 < -self.data.zero_tolerance {
            update.h_update_left += fwaves.1;
            update.hu_update_left += fwaves.1 * wavespeeds.1;
        } else {
            update.h_update_left += 0.5 * fwaves.1;
            update.hu_update_left += 0.5 * fwaves.1 * wavespeeds.1;
            update.h_update_right += 0.5 * fwaves.1;
            update.hu_update_right += 0.5 * fwaves.1 * wavespeeds.1;
        }

        update.max_wavespeed = f32::max(f32::abs(wavespeeds.0), f32::abs(wavespeeds.1));
        update
    }
}
