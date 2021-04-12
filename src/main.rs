#![allow(dead_code)]

use glacier::block::Block;

fn main() {
    println!("Hello, world!");
    //let a = gitter::Gitter{5,5};

    let mut a = Block::new(10, 10, 1.0, 1.0);
    a.compute_numerical_fluxes();
    println!("Computed stuff?");
}
