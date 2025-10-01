use std::borrow::Borrow;
use std::borrow::BorrowMut;
use std::ops::Add;
use std::ops::Mul;
use std::ops::Sub;

use num::Complex;

pub mod ring;

// TODO: Sparse polynomial is probably better off with other representation than a naive vector
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Polynomial<T> {
    pub coef: Vec<T>,
}

impl<T> Polynomial<T> {
    pub fn new(coef: Vec<T>) -> Self {
        Self { coef }
    }

    pub fn new_empty() -> Self {
        Self { coef: Vec::new() }
    }

    pub fn from_slice(coef: &[T]) -> Self
    where
        T: Clone,
    {
        let coef = coef.to_vec();

        Self { coef }
    }

    pub fn scale<V>(mut self, x: V) -> Self
    where
        T: Mul<V, Output = T> + Copy,
        V: Copy,
    {
        for c in &mut self.coef {
            *c = *c * x;
        }

        self
    }

    pub fn dot(&self, rhs: Self) -> Self
    where
        T: Mul<T, Output = T> + Copy,
    {
        if self.coef.len() != rhs.coef.len() {
            panic!("Mismatch length");
        }

        let mut product = Vec::new();
        for i in 0..self.coef.len() {
            product.push(self.coef[i] * rhs.coef[i]);
        }

        Polynomial::new(product)
    }

    pub fn coef_mul(&self, rhs: &Self) -> Self
    where
        T: Mul<T, Output = T> + Copy,
    {
        let mut c: Vec<_> = self
            .coef
            .iter()
            .zip(rhs.coef.iter())
            .map(|(a, b)| *a * *b)
            .collect();

        let l = std::cmp::min(self.coef.len(), rhs.coef.len());
        c.extend_from_slice(&self.coef[l..]);
        c.extend_from_slice(&rhs.coef[l..]);

        Self::new(c)
    }

    pub fn is_zero(&self) -> bool {
        self.coef.len() == 0
    }

    pub fn max_canon(&self) -> Option<&T>
    where
        T: Ord,
    {
        self.coef.iter().max()
    }
}

impl Polynomial<i64> {
    pub fn div_rem(&self, rhs: Polynomial<i64>) -> (Polynomial<i64>, Polynomial<i64>) {
        if self.is_zero() {
            return (Self::new_empty(), Self::new_empty());
        }

        let scale_val = rhs.coef[rhs.coef.len() - 1];
        let scaled_coefs: Vec<_> = rhs.coef[..rhs.coef.len() - 1].iter().cloned().collect();

        let divisor = Self::new(scaled_coefs).scale(1 / scale_val);
        if rhs.coef.len() > self.coef.len() {
            return (self.clone(), Self::new_empty());
        }

        let dim_result = self.coef.len() - rhs.coef.len();
        let mut rem = (*self).clone();
        let mut j = rem.coef.len() - 1;
        for i in (0..=dim_result).rev() {
            let result = Self::new(rem.coef[i..j].iter().cloned().collect())
                - divisor.clone().scale(rem.coef[j]);

            (&mut rem.coef[i..j]).copy_from_slice(&result.coef);
            j -= 1;
        }

        let res = Polynomial::from_slice(&rem.coef[j + 1..]);
        let rem = Polynomial::from_slice(&rem.coef[..j + 1]);

        (res.scale(1 / scale_val), rem)
    }

    pub fn translate(&mut self, q: i64) {
        self.coef.iter_mut().for_each(|e| *e = *e + q);
    }

    // Mod all the coefficient
    pub fn coef_mod(&mut self, q: i64) {
        self.coef.iter_mut().for_each(|e| *e = (*e).rem_euclid(q));
    }

    fn trim_zeroes(self) -> Polynomial<i64> {
        let trimmed: Vec<_> = self
            .coef
            .into_iter()
            .rev()
            .skip_while(|e| *e == 0)
            .collect();

        Self::new(trimmed)
    }
}

impl Polynomial<f64> {
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

        if self.coef.len() > other.coef.len() {
            sum.extend_from_slice(&self.coef[min_len..]);
        } else {
            sum.extend_from_slice(&other.coef[min_len..]);
        }

        Self { coef: sum }
    }
}

impl<T: Add<Output = T> + Copy> Add for &Polynomial<T> {
    type Output = Polynomial<T>;

    fn add(self, other: Self) -> Polynomial<T> {
        let min_len = std::cmp::min(self.coef.len(), other.coef.len());
        let mut sum = Vec::with_capacity(min_len);

        for i in 0..min_len {
            sum.push(self.coef[i] + other.coef[i]);
        }

        if self.coef.len() > other.coef.len() {
            sum.extend_from_slice(&self.coef[min_len..]);
        } else {
            sum.extend_from_slice(&other.coef[min_len..]);
        }

        Polynomial::new(sum)
    }
}

impl<T: Sub<Output = T> + Copy> Sub for Polynomial<T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let min_len = std::cmp::min(self.coef.len(), other.coef.len());
        let mut sum = Vec::with_capacity(min_len);

        for i in 0..min_len {
            sum.push(self.coef[i] - other.coef[i]);
        }

        if self.coef.len() > other.coef.len() {
            sum.extend_from_slice(&self.coef[min_len..]);
        } else {
            sum.extend_from_slice(&other.coef[min_len..]);
        }

        Self { coef: sum }
    }
}

impl<T: Mul<Output = T> + Copy + Add<Output = T> + Default> Mul<Polynomial<T>> for Polynomial<T> {
    type Output = Polynomial<T>;

    fn mul(self, other: Self) -> Self {
        let dim_result = self.coef.len() + other.coef.len();
        let mut coef = vec![T::default(); dim_result];
        for i in 0..self.coef.len() {
            for j in 0..other.coef.len() {
                coef[i + j] = coef[i + j] + self.coef[i] * other.coef[j];
            }
        }

        Self { coef }
    }
}

impl<'a, T: Mul<Output = T> + Copy + Add<Output = T> + Default> Mul<&'a Polynomial<T>>
    for &'a Polynomial<T>
{
    type Output = Polynomial<T>;

    fn mul(self, other: Self) -> Self::Output {
        let dim_result = self.coef.len() + other.coef.len();
        let mut coef = vec![T::default(); dim_result];
        for i in 0..self.coef.len() {
            for j in 0..other.coef.len() {
                coef[i + j] = coef[i + j] + self.coef[i] * other.coef[j];
            }
        }

        Polynomial::new(coef)
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

impl<T> core::ops::Index<usize> for Polynomial<T> {
    type Output = T;

    #[inline]
    fn index(&self, index: usize) -> &T {
        &self.coef[index]
    }
}

impl<T> Borrow<Vec<T>> for Polynomial<T> {
    fn borrow(&self) -> &Vec<T> {
        &self.coef
    }
}

impl<T> BorrowMut<Vec<T>> for Polynomial<T> {
    fn borrow_mut(&mut self) -> &mut Vec<T> {
        &mut self.coef
    }
}

pub fn gen_random_polynomial(length: usize, bound: usize) -> Polynomial<i64> {
    let mut coefs = vec![0; length];

    let bound = bound as i64;
    for i in 0..length {
        coefs[i] = rand::random_range(0..=bound);
    }

    Polynomial::from_slice(&coefs)
}

// /// Borrowed `Polynomial`.
// ///
// /// Similar to `Polynomial` that wraps a vector for its implementation, `PolyRef` wraps a slice.
// struct PolyRef(&[u8]);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_div_rem() {
        let p1 = Polynomial::from_slice(&[0, 0, 0, 5, 0, 0, 500]);
        let p2 = Polynomial::from_slice(&[1, 0, 0, 0, 1]);
        let (res, rem) = p1.div_rem(p2);

        let expected_res = Polynomial::from_slice(&[0, 0, 500]);
        assert_eq!(res, expected_res);

        let expected_rem = Polynomial::from_slice(&[0, 0, -500, 5]);
        assert_eq!(rem, expected_rem);
    }
}
