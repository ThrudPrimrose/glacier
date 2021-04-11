use super::gitter;

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
    boundaries: [BoundaryType; 4],
}

impl Block {
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
            boundaries: [BoundaryType::Wall; 4],
        }
    }
}
