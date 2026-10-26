use crate::poly::Polynomial;

/// Generate keys from R_2
pub fn uniform_R2(length: usize) -> Polynomial<i64> {
    let mut coefs = vec![0; length];

    for i in 0..length {
        coefs[i] = rand::random_range(-1..=1);
    }

    Polynomial::from_slice(&coefs)
}
