//! BFV Scheme
//!
//! Implementation reference: https://faculty.kfupm.edu.sa/coe/mfelemban/SEC595/References/Introduction%20to%20the%20BFV%20FHE%20Scheme.pdf
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
        let ek =
            EvalKey::generate_from_secret_key(&sk, self.params.Q, self.params.N, self.params.t);

        Session {
            params: self.params,
            sk,
            pk,
            ek,
        }
    }
}

pub struct Session {
    params: BGVParams,
    sk: SecretKey,
    pk: PublicKey,
    ek: EvalKey,
}

impl Session {
    pub fn encrypt(&self, m: Message) -> CipherText {
        let ring = self.ring_q();

        CipherText(ring.add(&m.0, &self.pk.0), self.pk.1.clone())
    }

    pub fn decrypt(&self, c: CipherText) -> Message {
        let ring = self.ring_q();

        let mut msg = c.0 + ring.mul(&c.1, &self.sk.0);
        msg.coef_mod(self.params.t as i64);

        Message(msg)
    }

    pub fn add_cipher(&self, lhs: CipherText, rhs: CipherText) -> CipherText {
        let ring = self.ring_q();

        CipherText(ring.add(&lhs.0, &rhs.0), ring.add(&lhs.1, &rhs.1))
    }

    pub fn mul_cipher(&self, lhs: CipherText, rhs: CipherText) -> CipherText {
        let ring = self.ring_q();

        let c0 = ring.mul(&lhs.0, &rhs.0);
        let c1 = ring.mul(&lhs.0, &rhs.1) + ring.mul(&lhs.1, &rhs.0);

        // TODO: so messy here
        // need something to do with i64 -> f64; maybe just use f64 instead?
        let c2 = ring.mul(&lhs.1, &rhs.1);

        let c0 = ring.add(&c0, &ring.mul(&self.ek.0, &c2));
        let c1 = ring.add(&c1, &ring.mul(&self.ek.1, &c2));

        let c0: Vec<i64> = c0
            .coef
            .into_iter()
            .map(|e| ((e * self.params.t as i64) as f64 / self.params.Q as f64).round() as i64)
            .collect();
        let c0 = Polynomial::new(c0);

        let c1: Vec<i64> = c1
            .coef
            .into_iter()
            .map(|e| ((e * self.params.t as i64) as f64 / self.params.Q as f64).round() as i64)
            .collect();
        let c1 = Polynomial::new(c1);

        CipherText(c0, c1)
    }

    fn ring_q(&self) -> ResidualQuotientRing {
        ResidualQuotientRing::cyclotomic(self.params.N, self.params.Q)
    }
}

#[derive(Debug, Clone)]
pub struct SecretKey(Polynomial<i64>);

#[derive(Debug, Clone)]
pub struct PublicKey(Polynomial<i64>, Polynomial<i64>);

impl PublicKey {
    fn generate_from_secret_key(s: &SecretKey, q: usize, N: usize, t: usize) -> Self {
        let a = gen_random_polynomial(N, q);
        let e = gen_random_polynomial(N, 1);

        let ring = ResidualQuotientRing::cyclotomic(N, q).ring;
        let pk_1 = ring.mul(&a, &s.0);

        Self(pk_1.scale(-1) + e.scale(t as i64), a)
    }
}

pub struct EvalKey(Polynomial<i64>, Polynomial<i64>);

impl EvalKey {
    fn generate_from_secret_key(s: &SecretKey, q: usize, N: usize, p: usize) -> Self {
        let a = gen_random_polynomial(N, q);
        let e = gen_random_polynomial(N, 1);
        let ring = ResidualQuotientRing::cyclotomic(N, q);

        let ek_1 = ring.add(&(ring.mul(&a, &s.0).scale(-1) + e), &ring.mul(&s.0, &s.0));

        Self(ek_1, a)
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
            t: 8,
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
    }

    #[test]
    fn test_addition() {
        let params = BGVParams {
            N: 4,
            t: 8,
            Q: 1024,
            L: 10,
            P: 1024,
        };
        let builder = SessionBuilder::new(params);
        let session = builder.generate_new();

        let input = Polynomial::from_slice(&[1, 2]);
        let msg = Message(input);
        let cipher_1 = session.encrypt(msg);
        let cipher_2 = cipher_1.clone();

        let c_sum = session.add_cipher(cipher_1, cipher_2);
        let decrypted = session.decrypt(c_sum);

        dbg!(&decrypted);

        assert!(false);
    }

    #[test]
    fn test_mul() {
        let params = BGVParams {
            N: 4,
            t: 8,
            Q: 1024,
            L: 10,
            P: 1024,
        };
        let builder = SessionBuilder::new(params);
        let session = builder.generate_new();

        let input = Polynomial::from_slice(&[1, 2]);
        let msg = Message(input);
        let cipher_1 = session.encrypt(msg);
        let cipher_2 = cipher_1.clone();

        let c_sum = session.mul_cipher(cipher_1, cipher_2);
        let decrypted = session.decrypt(c_sum);

        dbg!(&decrypted);

        assert!(false);
    }
}
