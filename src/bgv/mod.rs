use crate::ckks::codec::Codec;
use num::Complex;

use crate::poly::gen_random_polynomial;
use crate::poly::ring::*;
use crate::poly::Polynomial;

pub struct BGVParams {
    N: usize,
    t: usize,
    Q: usize,
    L: usize,
    P: usize,
}

pub struct SessionBuilder {
    params: BGVParams,
}

impl SessionBuilder {
    fn new(params: BGVParams) -> Self {
        Self { params }
    }

    pub fn generate_new(self) -> Session {
        let sk = SecretKey(gen_R2(self.params.N));
        let pk =
            PublicKey::generate_from_secret_key(&sk, self.params.Q, self.params.N, self.params.t);

        Session {
            params: self.params,
            sk,
            pk,
        }
    }
}

pub struct Session {
    params: BGVParams,
    sk: SecretKey,
    pk: PublicKey,
}

impl Session {
    pub fn encrypt(&self, m: Message) -> CipherText {
        let ring = self.ring();

        CipherText(ring.add(&m.0, &self.pk.0), self.pk.1.clone())
    }

    pub fn decrypt(&self, c: CipherText) -> Message {
        let ring = self.ring();

        Message(c.0 + ring.mul(&c.1, &self.sk.0))
    }

    pub fn add_cipher(&self, lhs: CipherText, rhs: CipherText) -> CipherText {
        let ring = self.ring();

        CipherText(ring.add(&lhs.0, &rhs.0), ring.add(&lhs.1, &rhs.1))
    }

    pub fn mul_cipher(&self, lhs: CipherText, rhs: CipherText) -> CipherText {
        todo!()
    }

    fn ring(&self) -> ResidualQuotientRing {
        ResidualQuotientRing::cyclotomic(self.params.N, self.params.t)
    }
}

#[derive(Debug, Clone)]
pub struct SecretKey(Polynomial<i64>);

#[derive(Debug, Clone)]
pub struct PublicKey(Polynomial<i64>, Polynomial<i64>);

impl PublicKey {
    fn generate_from_secret_key(s: &SecretKey, q: usize, N: usize, t: usize) -> Self {
        let a = gen_random_polynomial(N, q);
        let e = gen_random_polynomial(1, 1).scale(t as i64);
        let ring = ResidualQuotientRing::cyclotomic(N, t);
        let pk_1 = ring.mul(&a, &s.0);

        Self(pk_1.scale(-1) + e, a)
    }
}

#[derive(Debug, Clone)]
pub struct CipherText(Polynomial<i64>, Polynomial<i64>);

#[derive(Debug, Clone)]
pub struct Message(Polynomial<i64>);

/// Generate keys from R_2 with integer coefficients in {-1, 0, 1}
pub fn gen_R2(length: usize) -> Polynomial<i64> {
    let mut coefs = vec![0; length];

    for i in 0..length {
        coefs[i] = rand::random_range(-1..=1);
    }

    Polynomial::from_slice(&coefs)
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::cmp::PartialOrd;
    use std::ops::Sub;

    use num::Complex;

    #[test]
    fn test_encrypt_decrypt() {
        let params = BGVParams {
            N: 4,
            t: 1024,
            Q: 1024,
            L: 10,
            P: 1024,
        };
        let builder = SessionBuilder::new(params);
        let session = builder.generate_new();

        let input = Polynomial::from_slice(&[5, 4]);
        let msg = Message(input);
        let encrypted = session.encrypt(msg);
        let decrypted = session.decrypt(encrypted);

        assert!(false);

        // assert_complex_eq(decoded.0.into_iter(), expected.into_iter(), 0.001);
    }

    #[test]
    fn test_addition() {
        let params = BGVParams {
            N: 4,
            t: 1,
            Q: 1024,
            L: 10,
            P: 1024,
        };
        let builder = SessionBuilder::new(params);
        let session = builder.generate_new();

        let input = Polynomial::from_slice(&[5, 4]);
        let msg = Message(input);
        let cipher_1 = session.encrypt(msg);
        let cipher_2 = cipher_1.clone();

        let c_sum = session.add_cipher(cipher_1, cipher_2);
        let decrypted = session.decrypt(c_sum);

        dbg!(&decrypted);

        assert!(false);
    }

    fn assert_complex_eq<T>(
        x: impl ExactSizeIterator<Item = Complex<T>>,
        y: impl ExactSizeIterator<Item = Complex<T>>,
        epsilon: T,
    ) where
        T: Sub + PartialOrd + Clone + num::Num + num::Float,
        <T as Sub>::Output: PartialOrd<T>,
    {
        assert_eq!(x.len(), y.len(), "Object have different length");

        x.into_iter()
            .zip(y.into_iter())
            .for_each(|(a, b)| assert!((a - b).norm() < epsilon));
    }
}
