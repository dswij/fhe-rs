pub struct Polynomial<T> {
    coef: Vec<T>,
}

impl<T> Polynomial<T> {
    pub fn new(coef: Vec<T>) -> Self {
        Self { coef }
    }
}
