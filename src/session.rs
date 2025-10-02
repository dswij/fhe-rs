use crate::ckks::codec::Codec;
use num::Complex;

use crate::poly::gen_random_polynomial;
use crate::poly::ring::*;
use crate::poly::Polynomial;

// TODO: Bootstrapping operations?
// TODO: need to add a Bound
// TODO: need to cleanup the Polynomial: implement in f64 maybe instead of i64?
// TODO: add a borrowed polynomial
// TODO: use Pow instead of copying/referencing polynomial

pub struct CKKSParams {
    poly_degree: usize,
    q_0: usize,
    delta: usize,
    L: usize,
    P: usize,
}

pub struct SessionBuilder {
    params: CKKSParams,
}

impl SessionBuilder {
    fn new(params: CKKSParams) -> Self {
        Self { params }
    }

    pub fn from_existing_sk(self, sk: SecretKey) -> Session {
        let pk = PublicKey::generate_from_secret_key(&sk, self.params.q_0, self.params.poly_degree);
        let ek = EvalKey::generate_from_secret_key(
            &sk,
            self.params.q_0,
            self.params.poly_degree,
            self.params.P,
        );

        Session {
            params: self.params,
            sk,
            pk,
            ek,
        }
    }
}

pub struct Session {
    params: CKKSParams,
    sk: SecretKey,
    pk: PublicKey,
    ek: EvalKey,
}

impl Session {
    pub fn encode(&self, pt: PlainText) -> Message {
        let inner = self.codec().encode(pt.0);

        Message(inner)
    }

    pub fn decode(&self, m: Message) -> PlainText {
        let inner = self.codec().decode(m.0);

        PlainText(inner)
    }

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
        let ring = self.ring();
        let c0 = ring.mul(&lhs.0, &rhs.0);
        let c1 = ring.mul(&lhs.0, &rhs.1) + ring.mul(&lhs.1, &rhs.0);

        // TODO: so messy here
        // need something to do with i64 -> f64; maybe just use f64 instead?
        let c2 = ring.mul(&lhs.1, &rhs.1);
        let c2 = (
            Polynomial::new(
                ring.mul(&self.ek.0, &c2)
                    .coef
                    .into_iter()
                    .map(|e| ((e as f64) / self.params.P as f64).round() as i64)
                    .collect(),
            ),
            Polynomial::new(
                ring.mul(&self.ek.1, &c2)
                    .coef
                    .into_iter()
                    .map(|e| ((e as f64) / self.params.P as f64).round() as i64)
                    .collect(),
            ),
        );

        let c0: Vec<i64> = (c0 + c2.0)
            .coef
            .into_iter()
            .map(|e| ((e as f64) / self.params.delta as f64).round() as i64)
            .collect();
        let c0 = Polynomial::new(c0);

        let c1: Vec<i64> = (c1 + c2.1)
            .coef
            .into_iter()
            .map(|e| ((e as f64) / self.params.delta as f64).round() as i64)
            .collect();
        let c1 = Polynomial::new(c1);

        CipherText(c0, c1)
    }

    fn ring(&self) -> ResidualQuotientRing {
        ResidualQuotientRing::cyclotomic(self.params.poly_degree, self.params.q_0)
    }

    fn codec(&self) -> Codec {
        Codec::new(self.params.poly_degree, self.params.delta as f64)
    }
}

#[derive(Debug, Clone)]
pub struct SecretKey(Polynomial<i64>);

impl SecretKey {
    // TODO: fix this -- need some length cleanup
    pub fn generate_random() -> Self {
        Self(gen_random_polynomial(1, 1))
    }
}

#[derive(Debug, Clone)]
pub struct PublicKey(Polynomial<i64>, Polynomial<i64>);

impl PublicKey {
    fn generate_from_secret_key(s: &SecretKey, q: usize, N: usize) -> Self {
        let a = gen_random_polynomial(N, q);
        let ring = ResidualQuotientRing::cyclotomic(N, q);
        let pk_1 = ring.mul(&a, &s.0);

        Self(pk_1.scale(-1), a)
    }
}

#[derive(Debug, Clone)]
pub struct CipherText(Polynomial<i64>, Polynomial<i64>);

#[derive(Debug, Clone)]
pub struct Message(Polynomial<i64>);

#[derive(Debug, Clone)]
pub struct PlainText(Vec<Complex<f64>>);

pub struct EvalKey(Polynomial<i64>, Polynomial<i64>);

impl EvalKey {
    fn generate_from_secret_key(s: &SecretKey, q: usize, N: usize, p: usize) -> Self {
        let a = gen_random_polynomial(N, q);
        let ring = ResidualQuotientRing::cyclotomic(N, q);

        let ek_1 = ring.add(
            &ring.mul(&a, &s.0).scale(-1),
            &ring.mul(&s.0, &s.0).scale(p as i64),
        );

        Self(ek_1, a)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::cmp::PartialOrd;
    use std::ops::Sub;

    use num::Complex;

    #[test]
    fn test_encrypt_decrypt() {
        let params = CKKSParams {
            poly_degree: 4,
            q_0: 2048,
            delta: 64,
            L: 10,
            P: 256,
        };
        let builder = SessionBuilder::new(params);

        let sk = SecretKey(Polynomial::from_slice(&[1, 0, 0, 100]));
        let session = builder.from_existing_sk(sk);

        let input = vec![Complex::new(3.0, 4.0), Complex::new(2.0, -1.0)];
        let pt = PlainText(input);
        let msg = session.encode(pt);

        let encrypted = session.encrypt(msg);

        let decrypted = session.decrypt(encrypted);

        let decoded = session.decode(decrypted);
        let expected = vec![
            Complex::new(3.0082329989778316, 4.002601910021414),
            Complex::new(1.9917670010221684, -0.9973980899785868),
        ];

        assert_complex_eq(decoded.0.into_iter(), expected.into_iter(), 0.001);
    }

    #[test]
    fn test_addition() {
        let params = CKKSParams {
            poly_degree: 4,
            q_0: 1024,
            delta: 64,
            L: 10,
            P: 256,
        };
        let builder = SessionBuilder::new(params);

        let sk = SecretKey(Polynomial::from_slice(&[1, 0, 1, 0]));
        let session = builder.from_existing_sk(sk);

        let input = vec![Complex::new(3.0, 4.0), Complex::new(2.0, -1.0)];
        let pt = PlainText(input);

        let msg1 = session.encode(pt);
        let msg2 = msg1.clone();

        let c1 = session.encrypt(msg1);
        let c2 = session.encrypt(msg2);

        let c_sum = session.add_cipher(c1, c2);

        let decrypted = session.decrypt(c_sum);
        let decoded = session.decode(decrypted);

        let expected = vec![
            Complex::new(6.016465997955663, 8.005203820042828),
            Complex::new(3.983534002044337, -1.9947961799571736),
        ];

        assert_complex_eq(decoded.0.into_iter(), expected.into_iter(), 0.001);
    }

    #[test]
    fn test_mul() {
        let params = CKKSParams {
            poly_degree: 4,
            q_0: 1024,
            delta: 64 * 64,
            L: 10,
            P: 1,
        };
        let builder = SessionBuilder::new(params);

        let sk = SecretKey(Polynomial::from_slice(&[1, 0, 0, 0]));
        let session = builder.from_existing_sk(sk);

        let input = vec![Complex::new(3.0, 4.0), Complex::new(2.0, -1.0)];
        let pt = PlainText(input);

        let msg1 = session.encode(pt);
        let msg2 = msg1.clone();

        let c1 = session.encrypt(msg1);
        let c2 = session.encrypt(msg2);

        let c_mul = session.mul_cipher(c1, c2);

        let decrypted = session.decrypt(c_mul);
        let decoded = session.decode(decrypted);

        let expected = vec![
            Complex::new(-6.999465913858007, 24.00114644388103),
            Complex::new(2.999465913858012, -3.9998301186189744),
        ];

        assert_complex_eq(decoded.0.into_iter(), expected.into_iter(), 0.001);
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
