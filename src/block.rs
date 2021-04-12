use rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelIterator,
};

use super::gitter;
use crate::solvers::{fwave, Solver, Update};

#[derive(PartialEq, PartialOrd, Copy, Clone)]
enum BoundaryType {
    Outflow,
    Wall,
    Inflow,
    Connect,
    Passive,
}

#[derive(PartialEq, PartialOrd, Copy, Clone)]
enum BoundaryEdge {
    BndLeft,
    BndRight,
    BndBottom,
    BndTop,
}

pub struct Block {
    nx: usize,
    ny: usize,
    dx: f32,
    dy: f32,
    h: gitter::Gitter,
    hu: gitter::Gitter,
    hv: gitter::Gitter,
    b: gitter::Gitter,
    h_net_updates: gitter::Gitter,
    hu_net_updates: gitter::Gitter,
    hv_net_updates: gitter::Gitter,
    //boundaries: [BoundaryType; 4],
    solver: fwave::FWave,
    max_timestep: f32,
    offset_x: f32,
    offset_y: f32,
}

impl Block {
    #[inline]
    fn calculate_gitter_index(&self, i: usize) -> (usize, usize) {
        let y = i % (self.nx + 2);
        let x = i - (y * (self.nx + 2));
        (y, x)
    }

    #[inline]
    fn calculate_flat_index(&self, y: usize, x: usize) -> usize {
        y * (self.nx + 2) + x
    }

    pub fn new(l_nx: usize, l_ny: usize, l_dx: f32, l_dy: f32) -> Block {
        Block {
            nx: l_nx,
            ny: l_ny,
            dx: l_dx,
            dy: l_dy,
            h: gitter::Gitter::new(l_nx + 2, l_ny + 2),
            hu: gitter::Gitter::new(l_nx + 2, l_ny + 2),
            hv: gitter::Gitter::new(l_nx + 2, l_ny + 2),
            b: gitter::Gitter::new(l_nx + 2, l_ny + 2),
            h_net_updates: gitter::Gitter::new(l_nx + 2, l_ny + 2),
            hu_net_updates: gitter::Gitter::new(l_nx + 2, l_ny + 2),
            hv_net_updates: gitter::Gitter::new(l_nx + 2, l_ny + 2),
            //boundaries: [BoundaryType::Wall; 4],
            solver: fwave::FWave::default(),
            max_timestep: 0.0,
            offset_x: 0.0,
            offset_y: 0.0,
        }
    }

    pub fn compute_numerical_fluxes(&mut self) {
        let dx_inv: f32 = 1.0 / self.dx;
        let dy_inv: f32 = 1.0 / self.dy;

        //calculate updates in x direction
        let updates_x: Vec<Update> = (1..(self.ny + 1))
            .into_par_iter()
            .flat_map(|i| {
                (1..(self.nx + 2))
                    .into_par_iter()
                    .map(move |j| -> (usize, usize) { (i, j) })
            })
            .map(|(i, j)| -> Update {
                self.solver.compute_net_updates(
                    self.h[[i, j - 1]],
                    self.h[[i, j]],
                    self.hu[[i - 1, j - 1]],
                    self.hu[[i, j]],
                    self.b[[i, j - 1]],
                    self.b[[i, j]],
                )
            })
            .collect();

        //calculats updates in y direction
        let updates_y: Vec<Update> = (1..(self.ny + 2))
            .into_par_iter()
            .flat_map(|i| {
                (1..(self.nx + 1))
                    .into_par_iter()
                    .map(move |j| -> (usize, usize) { (i, j) })
            })
            .map(|(i, j)| -> Update {
                self.solver.compute_net_updates(
                    self.h[[i - 1, j]],
                    self.h[[i, j]],
                    self.hv[[i - 1, j]],
                    self.hv[[i, j]],
                    self.b[[i - 1, j]],
                    self.b[[i, j]],
                )
            })
            .collect();

        //get max wavespeed
        let id = || -> f32 { 0.0f32 };
        let op = |upd1: f32, upd2: f32| -> f32 { f32::max(upd1, upd2) };

        //extract reduced max_wavespeed from updates x
        let max_wavespeed_x: f32 = updates_x
            .par_iter()
            .cloned()
            .map(|updates| updates.max_wavespeed)
            .reduce(id, op);

        //extract reduced max_wavespeed from updates y
        let max_wavespeed_y: f32 = updates_y
            .par_iter()
            .cloned()
            .map(|updates| updates.max_wavespeed)
            .reduce(id, op);

        //get max wavespeed of both x and y directions
        let max_wavespeed = f32::max(max_wavespeed_x, max_wavespeed_y);

        //calculate max time step
        //updates x <- (1..(self.ny + 2)) ,  (1..(self.nx + 1))
        if max_wavespeed > 0.000001 {
            self.max_timestep = f32::min(self.dx / max_wavespeed, self.dy / max_wavespeed);
            self.max_timestep *= 0.4f32;
        } else {
            self.max_timestep = f32::MAX;
        }

        //can't access self within the loop so save the elements as constants before and access later
        let len = self.h.len();
        let nxp1 = self.nx + 1;
        let nxp2 = self.nx + 2;
        let nyp1 = self.ny + 1;
        let nyp2 = self.ny + 2;
        let dt = self.max_timestep;
        let lastlinebegin = len - nxp2;

        //update h from updates x
        //ranges for updates_x :(1..(self.ny + 1)),  (1..(self.nx + 2))
        /*
            -----------
            | |     | |
            | |     | |
            | |     | |
            -----------
            Ghost cells only receive updates from their right neighbor if first cell, and from their left update if they are the last cell of the line
            the updates are calculated only for between 1..nx+2 meaning that every line has nx+1 many Updates, need to map 2d coordinates of the gitter
            to the right offset in the flat updates array. This means from the 1d index that descrives 2d coordinates of the gitter
            one needs to remove the first ghost cell and the first and last line and update accordingly.

        */

        self.h.par_iter_mut().enumerate().for_each(|(i, el)| {
            //do not update ghost cells with wrong indices, first and last cell of everyrow
            //not in the first line (ghost cell), not in the last inge (ghost cells), not the first cell in the row
            if i > nxp2 && i < lastlinebegin {
                let offset = i % nxp2;
                let off = (i - (i / nxp2)) - nxp1;

                //first cell
                if offset == 0 {
                    *el = dt * dx_inv * updates_x[off].h_update_left;
                } else if offset == nxp2 - 1 {
                    //last cell
                    *el = dt * dx_inv * updates_x[off - 1].h_update_right;
                } else {
                    //not a ghost cell
                    *el -= dt
                        * (dx_inv * updates_x[off - 1].h_update_right
                            + dx_inv * updates_x[off].h_update_left);
                }
            }
        });

        // update h from updates y
        // ranges for y is different (1..(self.ny + 2)), (1..(self.nx + 1))
        // Same concept as the x direction but all the lines becomes blocks now
        // left -> down, right -> up in y direction
        self.h.par_iter_mut().enumerate().for_each(|(i, el)| {
            if i > nxp2 && i < len - 1 {
                let offset = i % nxp2;
                //for all the lines before substract 2, for the current line substract 1
                let line = i / nxp2;
                let off = (i - (2 * (line - 1))) - 1 - nxp2;

                if offset != 0 && offset != nxp1 {
                    if line == 1 {
                        *el = dt * dx_inv * updates_y[off].h_update_left;
                    } else if line == nyp1 {
                        *el = dt * dx_inv * updates_y[off - 1].h_update_right;
                    } else {
                        *el -= dt
                            * (dx_inv * updates_y[off - 1].h_update_right
                                + dx_inv * updates_y[off].h_update_left);
                    }
                }
            }
        });

        //update hu from updates x

        //update hv form updates y

        //set hu = 0 if h = 0
        //TODO: bounds
        /*
        let href = &self.h;
        self.hu.par_iter_mut().enumerate().for_each(|(i, e)| {
            let h = href[i];
            if  f32::le(&h,&0.1f32) {
                *e = 0.0f32;
            }
        });
        */

        //set hv = 0 if h = 0

        //set h = 0 if h < 0
    }
}
