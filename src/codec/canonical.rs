use num::Complex;
use std::f64::consts::PI;

use nalg::DMatrix;
use nalg::LU;

use crate::poly::Polynomial;

#[derive(Debug, Clone)]
pub struct CanonicalCodec {
    N: usize,
    delta: f64,
}

impl CanonicalCodec {
    pub fn new(N: usize, delta: f64) -> Self {
        Self { N, delta }
    }

    pub fn encode(&self, b: Vec<Complex<f64>>) -> Polynomial<i64> {
        let pi_inverse_b = self.pi_inverse(b);
        let encoded = self.sigma_inverse(pi_inverse_b);

        // Scale and round
        let coefs = encoded
            .into_iter()
            .map(|e| (e.re * self.delta).round() as i64)
            .collect();

        Polynomial::new(coefs)
    }

    pub fn decode(&self, p: Polynomial<i64>) -> Vec<Complex<f64>> {
        let rescaled = p
            .coef
            .into_iter()
            .map(|e| (e as f64) / self.delta)
            .collect();
        let rescaled = Polynomial::new(rescaled);

        let z = self.sigma(rescaled);

        self.pi(z)
    }

    fn sigma_inverse(&self, b: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
        let xi = Complex::from_polar(1.0, PI / (self.N as f64));
        let vandermonde = vandermonde(xi, self.N);
        let lu = LU::new(DMatrix::from_row_slice(self.N, self.N, &vandermonde));
        let coef = lu.solve(&DMatrix::from_vec(self.N, 1, b)).unwrap();

        coef.data.into()
    }

    fn sigma(&self, p: Polynomial<f64>) -> Vec<Complex<f64>> {
        let xi = Complex::from_polar(1.0, PI / self.N as f64);

        (0..self.N)
            .map(|i| {
                let root = xi.powi((2 * i + 1).try_into().unwrap());
                p.eval(root)
            })
            .collect()
    }

    fn pi(&self, z: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
        let n = (self.N / 2).try_into().unwrap();

        z.into_iter().take(n).collect()
    }

    fn pi_inverse(&self, z: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
        let mut z = z;
        let z_hat: Vec<Complex<f64>> = z.clone().into_iter().rev().map(|e| e.conj()).collect();
        z.extend(z_hat);

        z
    }
}

pub fn vandermonde(xi: Complex<f64>, N: usize) -> Vec<Complex<f64>> {
    let mut matrix = Vec::new();

    for i in 0..N {
        let root = xi.powi((2 * i + 1).try_into().unwrap());
        for j in 0..N {
            matrix.push(root.powi(j.try_into().unwrap()));
        }
    }

    matrix
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    #[test]
    fn test_vandermonde() {
        let xi = Complex::from_polar(1.0, 2.0 * PI / 8.0);
        let result = vandermonde(xi, 4);
        let expected: Vec<Vec<_>> = vec![
            vec![
                "1+0j",
                "0.7071067811865476+0.7071067811865475j",
                "2.220446049250313e-16+1j",
                "-0.7071067811865474+0.7071067811865477j",
            ],
            vec![
                "1+0j",
                "-0.7071067811865474+0.7071067811865477j",
                "-4.440892098500626e-16-1j",
                "0.707106781186548+0.707106781186547j",
            ],
            vec![
                "1+0j",
                "-0.7071067811865479-0.7071067811865471j",
                "1.1102230246251565e-15+1j",
                "0.7071067811865464-0.7071067811865487j",
            ],
            vec![
                "1+0j",
                "0.707106781186547-0.707106781186548j",
                "-1.3877787807814457e-15-1j",
                "-0.707106781186549-0.707106781186546j",
            ],
        ];
        let expected: Vec<Complex<f64>> = expected
            .into_iter()
            .flat_map(|e| {
                e.into_iter()
                    .filter_map(|s| Complex::<f64>::from_str(s).ok())
            })
            .collect();

        assert_eq!(result, expected);
    }

    #[test]
    fn test_encode_decode() {
        let input = vec![Complex::new(3.0, 4.0), Complex::new(2.0, -1.0)];
        let codec = CanonicalCodec::new(4, 64.0);
        let encoded = codec.encode(input.clone());
        let decoded = codec.decode(encoded);

        // TODO: equal comparion in floating should not be here
        let expected = vec![
            Complex::new(3.0082329989778316, 4.002601910021414),
            Complex::new(1.9917670010221684, -0.9973980899785868),
        ];

        assert!(decoded == expected);
    }
}
