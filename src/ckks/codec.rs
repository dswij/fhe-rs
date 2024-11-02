use num::Complex;
use std::f64::consts::PI;
use std::str::FromStr;

use nalg::DMatrix;
use nalg::LU;

use crate::poly::Polynomial;

pub struct Codec {
    M: i64,
    scale: f64,
    xi: Complex<f64>,
    basis: Vec<Complex<f64>>,
}

impl Codec {
    pub fn new(M: i64, scale: f64) -> Self {
        let xi = Complex::from_polar(1.0, 2.0 * PI / (M as f64));
        let vandermonde = vandermonde(xi, M);

        let n = (M / 2) as usize;
        let basis = DMatrix::from_vec(n, n, vandermonde).transpose().data.into();

        Self {
            M,
            scale,
            xi,
            basis,
        }
    }

    pub fn encode(&self, b: Vec<i64>) -> Polynomial<Complex<f64>> {
        self.sigma_inverse(b)
    }

    pub fn decode(&self, p: Polynomial<Complex<f64>>) -> Vec<Complex<f64>> {
        self.sigma(p)
    }

    fn sigma_inverse(&self, b: Vec<i64>) -> Polynomial<Complex<f64>> {
        let n: usize = (self.M / 2) as usize;
        let b: Vec<_> = b.into_iter().map(|e| Complex::new(e as f64, 0.0)).collect();
        let xi = Complex::from_polar(1.0, 2.0 * PI / (self.M as f64));
        let vandermonde = vandermonde(xi, self.M);
        let lu = LU::new(DMatrix::from_row_slice(n, n, &vandermonde));
        let coef = lu.solve(&DMatrix::from_vec(n, 1, b)).unwrap();

        Polynomial::new(coef.data.into())
    }

    fn sigma(&self, p: Polynomial<Complex<f64>>) -> Vec<Complex<f64>> {
        let xi = Complex::from_polar(1.0, 2.0 * PI / self.M as f64);

        (0..self.M / 2)
            .map(|i| {
                let root = xi.powi((2 * i + 1).try_into().unwrap());
                p.eval(root)
            })
            .collect()
    }

    fn pi(&self, z: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
        let n = (self.M / 4).try_into().unwrap();

        z.into_iter().take(n).collect()
    }

    fn pi_inverse(&self, z: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
        let mut z = z;
        let z_hat: Vec<Complex<f64>> = z.clone().into_iter().map(|e| e.conj()).collect();
        z.extend(z_hat);

        z
    }
}

pub fn vandermonde(xi: Complex<f64>, M: i64) -> Vec<Complex<f64>> {
    let N = M / 2;
    let mut matrix = Vec::new();

    for i in 0..N {
        let root = xi.powi((2 * i + 1).try_into().unwrap());
        for j in 0..N {
            matrix.push(root.powi(j.try_into().unwrap()));
        }
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
        let input = vec![1, 2, 3, 4];
        let codec = Codec::new(8, 1.0);
        let encoded = codec.encode(input.clone());
        let decoded = codec.decode(encoded);

        let re_decoded: Vec<_> = decoded.into_iter().map(|e| e.re.round() as i64).collect();

        assert!(re_decoded == input);
    }

    #[test]
    fn test_add() {
        let input_1 = vec![1, 2, 3, 4];
        let input_2 = vec![2, 3, 4, 5];
        let codec = Codec::new(8, 1.0);
        let sum = codec.encode(input_1) + codec.encode(input_2);
        let decoded = codec.decode(sum);
        let re_decoded: Vec<_> = decoded.into_iter().map(|e| e.re.round() as i64).collect();

        assert_eq!(re_decoded, [3, 5, 7, 9]);
    }

    #[test]
    fn test_mult() {
        let input_1 = vec![1, 2, 3, 4];
        let input_2 = vec![2, 3, 4, 5];
        let codec = Codec::new(8, 1.0);
        let sum = codec.encode(input_1) * codec.encode(input_2);
        let decoded = codec.decode(sum);
        let re_decoded: Vec<_> = decoded.into_iter().map(|e| e.re.round() as i64).collect();

        assert_eq!(re_decoded, [2, 6, 12, 20]);
    }

    // #[test]
    // fn test_pi_inverse() {
    //     let z = vec![0, 1];
    //
    //     let codec = Codec::new(8, 1.0);
    //
    //
    //     assert_eq!(codec.pi_inverse(z), [0, 1, 1, 0]);
    // }
}
