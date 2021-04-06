extern crate glacier;
pub use glacier::gitter;

#[cfg(test)]
mod tests {
    use glacier::gitter::Gitter;

    #[test]
    fn sanity_check() {
        assert_eq!(4, 2 + 2);
    }

    fn init_gitter() -> Gitter {
        let mut x = Gitter::new(10, 10);
        for i in 0..10 {
            for j in 0..10 {
                x[[i, j]] = (i * j) as f32;
            }
        }
        return x;
    }

    #[test]
    fn gitter_content() {
        let g = init_gitter();
        assert_eq!(g[[9, 9]], 81.0);
        assert_eq!(g[[8, 9]], 72.0);
        assert_eq!(g[[0, 9]], 0.0);
        assert_eq!(g[[5, 4]], 20.0);
    }

    #[test]
    fn get_row() {
        let g = init_gitter();
        let gr = g.get_row(1);
        assert_eq!(gr[0], 0.0);
        assert_eq!(gr[5], 5.0);
        assert_eq!(gr[1], 1.0);
    }
}
