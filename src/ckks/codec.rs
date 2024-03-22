use num::Complex;
use std::f64::consts::PI;
use std::str::FromStr;

use crate::linarg::try_solve_complex;
use crate::poly::Polynomial;

pub struct Codec {
    M: i64,
}

impl Codec {
    pub fn new(M: i64) -> Self {
        Self { M }
    }

    // TODO: implement scaler etc
    pub fn encode(&self, b: Vec<i64>) -> Polynomial<Complex<f64>> {
        let b = b
            .into_iter()
            .map(|e| vec![Complex::new(e as f64, 0.0)])
            .collect();
        let xi = Complex::from_polar(1.0, 2.0 * PI / (self.M as f64));
        let vandermonde = vandermonde(xi, self.M);

        let coef = try_solve_complex(vandermonde, b)
            .unwrap()
            .into_iter()
            .flatten()
            .collect();

        Polynomial::new(coef)
    }

    pub fn decode(&self, p: Polynomial<Complex<f64>>) -> Vec<Complex<f64>> {
        let xi = Complex::from_polar(1.0, 2.0 * PI / self.M as f64);

        (0..self.M / 2)
            .map(|i| {
                let root = xi.powi((2 * i + 1).try_into().unwrap());
                p.eval(root)
            })
            .collect()
    }
}

pub fn vandermonde(xi: Complex<f64>, M: i64) -> Vec<Vec<Complex<f64>>> {
    let N = M / 2;
    let mut matrix = Vec::new();

    for i in 0..N {
        let root = xi.powi((2 * i + 1).try_into().unwrap());
        let mut row = Vec::new();

        for j in 0..N {
            row.push(root.powi(j.try_into().unwrap()));
        }

        matrix.push(row);
    }

    matrix
}

mod test {
    use super::*;

    #[test]
    fn test_vandermonde() {
        let xi = Complex::from_polar(1.0, 2.0 * PI / 8.0);
        let result = vandermonde(xi, 8);
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
        let expected: Vec<Vec<Complex<f64>>> = expected
            .into_iter()
            .map(|e| {
                e.into_iter()
                    .filter_map(|s| Complex::<f64>::from_str(s).ok())
                    .collect()
            })
            .collect();

        assert_eq!(result, expected);
    }

    // TODO: encode/decode test
    //
    #[test]
    fn test_encode() {
        let input = vec![1, 2, 3, 4];
        let codec = Codec::new(8);
        let encoded = codec.encode(input);

        dbg!(&encoded);

        let decoded = codec.decode(encoded);
        dbg!(&decoded);
        assert!(false);
    }

    #[test]
    fn test_decode() {}
}
