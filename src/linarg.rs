use lapacke::{dgesv, zgesv, Layout};
use num::Complex;
use num::FromPrimitive;

// TODO: use macros to make it generic i guess
pub fn map_domain(
    x: Vec<Complex<f64>>,
    old: Vec<Complex<f64>>,
    new: Vec<Complex<f64>>,
) -> Vec<Complex<f64>>
where
{
    let offset = (old[1] * new[0] - old[0] * new[1]) / (old[1] - old[0]);
    let scl = (new[1] - new[0]) / (old[1] - old[0]);

    x.into_iter().map(|e| offset + scl * e).collect()
}

// TODO: make this generic instead
pub fn try_solve_complex(
    A: Vec<Vec<Complex<f64>>>,
    B: Vec<Vec<Complex<f64>>>,
) -> Result<Vec<Vec<Complex<f64>>>, String> {
    let n: i32 = A[0].len().try_into().unwrap();
    let nrhs: i32 = B[0].len().try_into().unwrap();

    let mut a: Vec<_> = A.into_iter().flatten().collect();
    let mut b: Vec<_> = B.into_iter().flatten().collect();
    let mut x = vec![1; n as usize];

    dbg!(&a);
    dbg!(&b);

    let info: i32;
    unsafe {
        info = zgesv(Layout::RowMajor, n, nrhs, &mut a, n, &mut x, &mut b, nrhs);
    }

    if info == 0 {
        let nrhs: usize = nrhs as usize;
        let n: usize = n as usize;
        let mut result = Vec::with_capacity(n);
        for i in 0..n {
            let mut row = Vec::with_capacity(nrhs);
            for j in 0..nrhs {
                row.push(b[i as usize * nrhs + j as usize]);
            }

            result.push(row);
        }

        dbg!(&result);

        Ok(result)
    } else if info < 0 {
        Err(format!("Argument {} is invalid", info))
    } else {
        Err(format!("U({}, {}) is zero", info, info))
    }
}

pub fn try_solve(A: Vec<Vec<f64>>, B: Vec<Vec<f64>>) -> Result<Vec<Vec<f64>>, String> {
    let n: i32 = A[0].len().try_into().unwrap();
    let nrhs: i32 = B[0].len().try_into().unwrap();

    let mut a: Vec<_> = A.into_iter().flatten().collect();
    let mut b: Vec<_> = B.into_iter().flatten().collect();
    let mut x = vec![0; n as usize];

    unsafe {
        let info = dgesv(Layout::RowMajor, n, nrhs, &mut a, n, &mut x, &mut b, nrhs);
        if info == 0 {
            let nrhs: usize = nrhs as usize;
            let n: usize = n as usize;
            let mut result = Vec::with_capacity(n);
            for i in 0..n {
                let mut row = Vec::with_capacity(nrhs);
                for j in 0..nrhs {
                    row.push(b[i as usize * nrhs + j as usize]);
                }

                result.push(row);
            }

            Ok(result)
        } else if info < 0 {
            Err(format!("Argument {} is invalid", info))
        } else {
            Err(format!("U({}, {}) is zero", info, info))
        }
    }
}

mod test {
    use super::*;

    #[test]
    fn test_dgesv() {
        let A = vec![vec![1.0, 2.0], vec![3.0, 5.0]];
        let B = vec![vec![1.0], vec![2.0]];

        let result = try_solve(A, B).unwrap();
        let expected: Vec<f64> = vec![-1.0, 1.0];

        // TODO: use something better for float comparison
        result
            .into_iter()
            .flatten()
            .zip(expected.into_iter())
            .map(|(a, b)| {
                assert!((a - b).abs() < 1e-10);
            });
    }

    #[test]
    fn test_zgesv() {
        let A = vec![vec![1.0, 2.0], vec![3.0, 5.0]];
        let B = vec![vec![1.0], vec![2.0]];

        // TODO: monstrosity
        let A = A
            .into_iter()
            .map(|e| {
                e.into_iter()
                    .map(|e| Complex::<f64>::from_f64(e).unwrap())
                    .collect()
            })
            .collect();
        let B = B
            .into_iter()
            .map(|e| {
                e.into_iter()
                    .map(|e| Complex::<f64>::from_f64(e).unwrap())
                    .collect()
            })
            .collect();

        let result = try_solve_complex(A, B).unwrap();
        let expected: Vec<f64> = vec![-1.0, 1.0];

        result
            .into_iter()
            .flatten()
            .zip(expected.into_iter())
            .map(|(a, b)| {
                assert!((a - b).re.abs() < 1e-10);
            });
    }
}
