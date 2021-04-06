extern crate glacier;
pub use glacier::gitter;

#[cfg(test)]
mod tests {
    #[test]
    fn sanity_check() {
        assert_eq!(4, 2 + 2);
    }
}
