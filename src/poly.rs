use std::ops::Add;
use std::ops::Mul;

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
                acc + coef * x.powu(i.try_into().unwrap())
            })
    }
}

impl<T: Add<Output = T> + Copy> Add for Polynomial<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let min_len = std::cmp::min(self.coef.len(), other.coef.len());
        let mut sum = Vec::with_capacity(min_len);

        for i in 0..min_len {
            sum.push(self.coef[i] + other.coef[i]);
        }

        Self { coef: sum }
    }
}

impl<T: Mul<Output = T> + Copy> Mul for Polynomial<T> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let min_len = std::cmp::min(self.coef.len(), other.coef.len());
        let mut coef = Vec::with_capacity(min_len);

        for i in 0..min_len {
            coef.push(self.coef[i] * other.coef[i]);
        }

        Self { coef }
    }
}

impl<T: std::fmt::Display> std::fmt::Display for Polynomial<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let coefs_with_symbol = self
            .coef
            .iter()
            .enumerate()
            .map(|(i, coef)| format!("{}x^{}", coef, i))
            .collect::<Vec<String>>()
            .join(" + ");

        write!(f, "{}", coefs_with_symbol)
    }
}
