//! Encode a decimal integer input to a n-base decomposition.
//!
//! Unused bits are set to 0.
//!
//! For example, given a decimal input `19`, the decomposition encoding to binary with degree 8 is
//!
//! 19 (dec) = 0 x^8 + 0 x^7, + 0 x^6 + 0 x^5 + 1 x^4 + 0 x^3 + 0 x^2 + 1 x + 1

use crate::Polynomial;

pub struct IntegerCodec {
    base: i64,
    degree: usize,
}

impl IntegerCodec {
    pub fn new(base: i64, degree: usize) -> Self {
        Self { base, degree }
    }

    pub fn encode(&self, input: i64) -> Polynomial<i64> {
        let mut r = Vec::with_capacity(self.degree + 1);
        let mut quotient = input;
        let mut rem;
        while quotient > 0 {
            (quotient, rem) = (
                quotient.div_euclid(self.base),
                quotient.rem_euclid(self.base),
            );

            r.push(rem);
        }

        r.resize(r.capacity(), 0);

        Polynomial::new(r)
    }

    // Decoding means evaluating at `x = degree`.
    pub fn decode(&self, decoded: Polynomial<i64>) -> i64 {
        decoded
            .coef
            .into_iter()
            .enumerate()
            .fold(0, |acc, (i, val)| {
                dbg!(i as u32);
                dbg!(self.base.pow(i as u32));

                acc + val * self.base.pow(i as u32) as i64
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_integer() {
        let codec = IntegerCodec::new(2, 8);
        let encoded = codec.encode(19);
        let expected = Polynomial::from_slice(&[1, 1, 0, 0, 1, 0, 0, 0, 0]);

        assert_eq!(encoded, expected);
    }

    #[test]
    fn test_decode_integer() {
        let codec = IntegerCodec::new(2, 8);
        let encoded = Polynomial::from_slice(&[1, 1, 0, 0, 1, 0, 0, 0, 0]);
        let decoded = codec.decode(encoded);

        assert_eq!(decoded, 19);
    }
}
