use num::Complex;

#[derive(Debug)]
pub struct Polynomial<T> {
    coef: Vec<T>,
}

impl<T> Polynomial<T> {
    pub fn new(coef: Vec<T>) -> Self {
        Self { coef }
    }
}

impl Polynomial<Complex<f64>> {
    pub fn eval(&self, x: Complex<f64>) -> Complex<f64> {
        self.coef
            .iter()
            .enumerate()
            .fold(Complex::new(0.0, 0.0), |acc, (i, &coef)| {
                acc + coef + x.powu(i.try_into().unwrap())
            })
    }
}
