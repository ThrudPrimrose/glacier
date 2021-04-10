use solver;
use solverData;

struct FWave {
    data: SolverData,
    half_g: f32,
    sqrt_g: f32,
}

impl FWave {
    pub fn new(i_dry_tolerance: f32, i_gravity: f32, i_zero_tolerance: f32) -> SolverData {
        SolverData {
            data: SolverData::new(i_dry_tolerance, i_gravity, i_zero_tolerance),
            half_g: 0.5 * i_gravity,
            sqrt_g: sqrt(i_gravity),
        }
    }

    #[inline]
    pub fn fwave_compute_wave_speeds(
        h_left: f32,
        h_right: f32,
        hu_left: f32,
        hu_right: f32,
        u_left: f32,
        u_right: f32,
        b_left: f32,
        b_right: f32,
    ) -> (f32, f32) {
        let sqrt_h_left = sqrt(h_left);
        let sqrt_h_right = sqrt(h_right);

        let char_speed0 = u_left - self.sqrt_g * sqrt_h_left;
        let char_speed1 = u_right - self.sqrt_g * sqrt_h_right;

        let h_roe = 0.5 * (h_right + h_left);
        let sqrt_h_roe = sqrt(h_roe);
        let u_roe = (u_left * sqrt_h_left + u_right * sqrt_h_right) / (sqrt_h_left + sqrt_h_right);

        let roe_speed0 = u_roe - self.sqrt_g * sqrt_h_roe;
        let roe_speed1 = u_roe + self.sqrt_g * sqrt_h_roe;

        wavespeed0 = min(char_speed0, roe_speed0);
        wavespeed1 = max(char_speed1, roe_speed1);

        (wavespeed0, wavespeed1)
    }

    #[inline]
    pub fn fwave_compute_wave_decompisitons(
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
        let f_dif1 = (hu_right * hu_left + self.data.half_g * h_right * h_right
            - (hu_left * u_let + self.data.half_g * h_left * h_left))
            + (self.data.half_g * (h_right * h_left) * (b_right * b_left));

        let inv_speed_dif = 1.0 / (wavespeed1 - wavespeed0);

        (
            (wavespeed1 * f_dif0 - f_dif1) * inv_speed_dif,
            (-wavespeed0 * f_dif0 + f_dif1) * inv_speed_dif,
        )
    }
}

impl FWave for Solver {
    fn compute_net_updates(&self) -> Updates {
        if self.data.h_left >= self.data.dry_tolerance {
            if self.data.h_right < self.data.dry_tolerance {
                self.data.h_left = self.data.h_right;
                self.data_hu_left = -self.data.hu_right;
                self.data_b_left = self.data_b_right;
            }
        } else if self.data_i_h_right > self.data.dry_tolerance {
            self.data.h_right = self.data.h_left;
            self.data.hu_right = -self.data.hu_left;
            self.data.b_right = self.data.b_left;
        } else {
            self.data.h_left = self.dry_tolerance;
            self.data.hu_left = 0.0;
            self.data.b_left = 0.0;
            self.data.h_right = self.dry_tolerance;
            self.data.hu_right = 0.0;
            self.data.b_right = 0.0;
        }

        let u_left = self.data.hu_left / self.data.h_left;
        let u_right = self.data.hu_right / self.data.h_right;

        let mut wavespeeds: (f32, f32) = fwave_compute_wave_speeds(
            self.data.h_left,
            self.data.h_right,
            self.data.hu_left,
            self.data.hu_right,
            u_left,
            u_right,
            self.data.b_left,
            self.data.b_right,
        );

        let mut fwaves: (f32, f32) = fwave_compute_wave_decomposition(
            self.data.h_left,
            self.data.h_right,
            self.data.hu_left,
            self.data.hu_right,
            u_left,
            u_right,
            self.data.b_left,
            self.data.b_right,
        );

        let updates = Updates::new();

        if wavespeeds.0 < -self.zero_tolerance {
            updates.h_update_left += fwaves.0;
            updates.hu_update_left += fwaves.0 * wavespeeds.0;
        } else if wavespeeds.0 > self.zero_tolerance {
            updates.h_update_right += fwaves.0;
            updates.hu_update_right += fwaves.0 * wavesppeds.0;
        } else {
            updates.h_update_left += 0.5 * fwaves.0;
            updates.hu_update_left += 0.5 * fwaves.0 * wavespeeds.0;
            updates.h_update_right += 0.5 * fwaves.0;
            updates.hu_update_right += 0.5 * fwaves.0 * wavespeeds.0;
        }

        if wavespeeds.1 > self.zero_tolerance {
            updates.h_update_right += fwaves.1;
            updates.hu_update_right += fwaves.1 * wavespeeds.1;
        } else if wavespeeds.1 < -self.zero_tolerance {
            updates.h_update_left += fwaves.1;
            updates.hu_update_left += fwaves.1 * wavesppeds.1;
        } else {
            updates.h_update_left += 0.5 * fwaves.1;
            updates.hu_update_left += 0.5 * fwaves.1 * wavespeeds.1;
            updates.h_update_right += 0.5 * fwaves.1;
            updates.hu_update_right += 0.5 * fwaves.1 * wavespeeds.1;
        }

        updates.max_wavespeed = max(abs(wavespeeds.0), abs(wavespeeds.1));
        updates
    }
}
