use std::{
    ops::{Index, IndexMut},
    slice::Iter,
};

use rayon::{iter::IntoParallelRefMutIterator, slice::IterMut};

pub struct Gitter {
    x: usize,         //size in x
    y: usize,         //size in y
    gitter: Vec<f32>, //use vector because array always require compile time constants
}

impl Gitter {
    pub fn new(x: usize, y: usize) -> Gitter {
        Gitter {
            x: x,
            y: y,
            gitter: vec![0.0; x * y],
        }
    }

    pub fn get_rows(&self) -> usize {
        self.y
    }

    pub fn get_cols(&self) -> usize {
        self.x
    }

    pub fn get_row(&self, row: usize) -> &[f32] {
        &self.gitter[(row * self.x)..(row * self.x + self.x)]
    }

    pub fn iter(&mut self) -> Iter<f32> {
        self.gitter.iter()
    }

    pub fn par_iter_mut(&mut self) -> IterMut<f32> {
        self.gitter.par_iter_mut()
    }

    //non-cont slice can only be a copy, leave it out for now
    //pub fn get_col(&self, row : usize) -> [f32] {}
}

impl Index<[usize; 2]> for Gitter {
    type Output = f32;

    fn index(&self, idx: [usize; 2]) -> &f32 {
        &self.gitter[idx[1] * self.x + idx[0]]
    }
}

impl IndexMut<[usize; 2]> for Gitter {
    fn index_mut(&mut self, idx: [usize; 2]) -> &mut f32 {
        &mut self.gitter[idx[1] * self.x + idx[0]]
    }
}

impl Index<usize> for Gitter {
    type Output = f32;

    fn index(&self, idx: usize) -> &f32 {
        &self.gitter[idx]
    }
}

impl IndexMut<usize> for Gitter {
    fn index_mut(&mut self, idx: usize) -> &mut f32 {
        &mut self.gitter[idx]
    }
}
