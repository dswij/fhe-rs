pub struct BFVParams {
    N: usize,
    t: usize,
    Q: usize,
    P: usize,
}

impl BFVParams {
    fn delta(&self) -> i64 {
        (self.Q / self.t) as i64
    }
}
