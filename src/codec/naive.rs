//! Naive Codec encodes a plaintext to the constant term of a polynomial.
//!
//! i.e. for an encoding to a polynomial of degree `N`,  the input `c` is encoded as
//! M = c + 0 x + 0 x^2 + ... + 0 x^{n-1}

use crate::Polynomial;

pub struct NaiveCodec {
    degree: usize,
}

impl NaiveCodec {
    pub fn new(degree: usize) -> Self {
        Self { degree }
    }

    pub fn encode(&self, input: i64) -> Polynomial<i64> {
        let mut coefs = vec![0; self.degree + 1];
        coefs[0] = input;

        Polynomial::new(coefs)
    }

    pub fn decode(&self, encoded: Polynomial<i64>) -> i64 {
        encoded[0]
    }
}

// This is so simple we might not need this.
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode() {
        let codec = NaiveCodec::new(8);
        let encoded = codec.encode(19);
        let expected = Polynomial::from_slice(&[19, 0, 0, 0, 0, 0, 0, 0, 0]);

        assert_eq!(encoded, expected);
    }

    #[test]
    fn test_decode() {
        let codec = NaiveCodec::new(8);
        let encoded = Polynomial::from_slice(&[19, 0, 0, 0, 0, 0, 0, 0, 0]);
        let decoded = codec.decode(encoded);

        assert_eq!(decoded, 19);
    }
}
