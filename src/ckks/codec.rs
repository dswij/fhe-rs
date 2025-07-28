use num::Complex;
use std::f64::consts::PI;
use std::str::FromStr;

use nalg::DMatrix;
use nalg::DVector;
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
        let b: Vec<_> = b.into_iter().map(|e| Complex::new(e as f64, 0.0)).collect();
        let pi_inverse_b = self.pi_inverse(b);
        let scaled: Vec<Complex<f64>> = pi_inverse_b.into_iter().map(|e| e * self.scale).collect();
        let proj = self.project_to_lattice(scaled);
        let encoded = self.sigma_inverse(proj);

        Polynomial::new(
            encoded
                .into_iter()
                .map(|e| Complex::new(e.re.round(), 0.0))
                .collect(),
        )
    }

    pub fn decode(&self, mut p: Polynomial<Complex<f64>>) -> Vec<Complex<f64>> {
        let rescaled = p.scale(1.0 / self.scale);
        let z = self.sigma(p);

        self.pi(z)
    }

    fn sigma_inverse(&self, b: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
        let n: usize = (self.M / 2) as usize;
        let xi = Complex::from_polar(1.0, 2.0 * PI / (self.M as f64));
        let vandermonde = vandermonde(xi, self.M);
        let lu = LU::new(DMatrix::from_row_slice(n, n, &vandermonde));
        let coef = lu.solve(&DMatrix::from_vec(n, 1, b)).unwrap();

        coef.data.into()
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
        let z_hat: Vec<Complex<f64>> = z.clone().into_iter().rev().map(|e| e.conj()).collect();
        z.extend(z_hat);

        z
    }

    fn compute_on_basis(&self, z: Vec<Complex<f64>>) -> Vec<f64> {
        let size = (self.M / 2) as usize;
        let basis_matrix = DMatrix::from_row_slice(size, size, &self.basis);
        let mut output = Vec::with_capacity(self.basis.len());
        let z = DVector::from_vec(z);

        for base in basis_matrix.row_iter() {
            let e = z.conjugate_transpose().dot(&base) / base.conjugate().dot(&base);
            output.push(e.re);
        }

        output
    }

    fn project_to_lattice(&self, z: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
        dbg!(&z);
        let size = (self.M / 2) as usize;
        assert_eq!(z.len(), size);

        let coords: Vec<_> = self
            .compute_on_basis(z)
            .into_iter()
            .map(|e| Complex::new(e.round(), 0.0))
            .collect();
        let coords_vector = DVector::from_vec(coords);

        let basis_matrix = DMatrix::from_row_slice(size, size, &self.basis);

        basis_matrix.tr_mul(&coords_vector).data.into()
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
        let input = vec![3, 2];
        let codec = Codec::new(8, 64.0);
        let encoded = codec.encode(input.clone());
        let decoded = codec.decode(encoded);

        let re_decoded: Vec<_> = decoded.into_iter().map(|e| e.re.round() as i64).collect();

        assert!(re_decoded == input);
    }

    #[test]
    fn test_add() {
        let input_1 = vec![3, 4];
        let input_2 = vec![4, 5];
        let codec = Codec::new(8, 64.0);
        let sum = codec.encode(input_1) + codec.encode(input_2);
        let decoded = codec.decode(sum);
        let re_decoded: Vec<_> = decoded.into_iter().map(|e| e.re.round() as i64).collect();

        assert_eq!(re_decoded, [7, 9]);
    }

    #[test]
    fn test_compute_basis() {
        let input = vec!["192+256j", "128-64j", "128+64j", "192-256j"]
            .into_iter()
            .filter_map(|s| Complex::<f64>::from_str(s).ok())
            .collect();
        let codec = Codec::new(8, 64.0);
        let result = codec.compute_on_basis(input);

        // TODO: equal comparison in floating number, really?
        assert_eq!(
            result,
            [
                160.0,
                90.5096679918781,
                159.99999999999997,
                45.25483399593886,
            ]
        );
    }

    // TODO: not working?
    #[test]
    fn test_mult() {
        let input_1 = vec![1, 2, 3, 4];
        let input_2 = vec![2, 3, 4, 5];
        let codec = Codec::new(8, 64.0);
        let sum = codec.encode(input_1) * codec.encode(input_2);
        let decoded = codec.decode(sum);
        let re_decoded: Vec<_> = decoded.into_iter().map(|e| e.re.round() as i64).collect();

        assert_eq!(re_decoded, [2, 6, 12, 20]);
    }

    // TODO: write out unit test for everything.
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
