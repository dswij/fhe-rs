use crate::ckks::codec::Codec;
use num::Complex;

use crate::poly::gen_random_polynomial;
use crate::poly::ring::*;
use crate::poly::Polynomial;

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
    fn encode(&self, pt: PlainText) -> Message {
        let codec = Codec::new(self.params.poly_degree, self.params.delta as f64);
        let inner = codec.encode(pt.0);

        Message(inner)
    }

    fn decode(&self, m: Message) -> PlainText {
        let codec = Codec::new(self.params.poly_degree, self.params.delta as f64);
        let inner = codec.decode(m.0);

        PlainText(inner)
    }

    fn encrypt(&self, m: Message) -> CipherText {
        let ring = ResidualQuotientRing::cyclotomic(self.params.poly_degree, self.params.q_0);

        CipherText(ring.add(&m.0, &self.pk.0), self.pk.1.clone())
    }

    fn decrypt(&self, c: CipherText) -> Message {
        // TODO: probably just store the ring too
        let ring = ResidualQuotientRing::cyclotomic(self.params.poly_degree, self.params.q_0);

        Message(c.0 + ring.mul(&c.1, &self.sk.0))
    }

    fn add_cipher(&self, lhs: CipherText, rhs: CipherText) -> CipherText {
        let ring = ResidualQuotientRing::cyclotomic(self.params.poly_degree, self.params.q_0);

        CipherText(ring.add(&lhs.0, &rhs.0), ring.add(&lhs.1, &rhs.1))
    }

    fn mul_cipher(&self, lhs: CipherText, rhs: CipherText) -> CipherText {
        println!("lhs: {}", lhs.0);
        println!("rhs: {}", rhs.0);
        println!("lhs * rhs: {}", &lhs.0 * &rhs.0);
        let ring = ResidualQuotientRing::cyclotomic(self.params.poly_degree, self.params.q_0);
        println!("quotient: {}", ring.ring.quotient);
        let c0 = ring.mul(&lhs.0, &rhs.0);
        println!("c__0: {}", c0);
        let c1 = ring.mul(&lhs.0, &rhs.1) + ring.mul(&lhs.1, &rhs.0);
        println!("c__1: {}", c1);

        // TODO: so messy here
        // need something to do with i64 -> f64; maybe just use f64 instead?
        let c2 = ring.mul(&lhs.1, &rhs.1);
        println!("c__2: {}", c2);
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
        println!("c_mul_0: {}", c0);

        let c1: Vec<i64> = (c1 + c2.1)
            .coef
            .into_iter()
            .map(|e| ((e as f64) / self.params.delta as f64).round() as i64)
            .collect();
        let c1 = Polynomial::new(c1);
        println!("c_mul_1: {}", c1);

        CipherText(c0, c1)
    }

    fn mul_cipher_and_decrypt(&self, lhs: CipherText, rhs: CipherText) -> Message {
        let ring = ResidualQuotientRing::cyclotomic(self.params.poly_degree, self.params.q_0);
        let c0 = ring.mul(&lhs.0, &rhs.0);
        let c1 = ring.mul(&lhs.0, &rhs.1) + ring.mul(&lhs.1, &rhs.0);

        // TODO: so messy here
        // need something to do with i64 -> f64; maybe just use f64 instead?
        let c2 = ring.mul(&lhs.1, &rhs.1);

        let result =
            c0 + ring.mul(&c1, &self.sk.0) + ring.mul(&c2, &ring.mul(&self.sk.0, &self.sk.0));

        Message(Polynomial::new(
            result
                .coef
                .into_iter()
                .map(|e| ((e as f64) / self.params.delta as f64).round() as i64)
                .collect(),
        ))
    }
}

#[derive(Debug, Clone)]
pub struct SecretKey(Polynomial<i64>);

impl SecretKey {
    fn generate_random() -> Self {
        Self(gen_random_polynomial(1, 1))
    }
}

#[derive(Debug, Clone)]
pub struct PublicKey(Polynomial<i64>, Polynomial<i64>);

impl PublicKey {
    fn generate_from_secret_key(s: &SecretKey, q: usize, N: usize) -> Self {
        let a = Polynomial::from_slice(&[1, 0, 0, 0]);
        let ring = ResidualQuotientRing::cyclotomic(N, q);
        // let pk_1 = ring.add(&ring.mul(&a, &s.0), &e);
        let pk_1 = ring.mul(&a, &s.0);

        Self(pk_1.scale(-1), a)
    }
}

#[derive(Debug, Clone)]
struct CipherText(Polynomial<i64>, Polynomial<i64>);

#[derive(Debug, Clone)]
struct Message(Polynomial<i64>);

#[derive(Debug, Clone)]
struct PlainText(Vec<Complex<f64>>);

struct EvalKey(Polynomial<i64>, Polynomial<i64>);

impl EvalKey {
    fn generate_from_secret_key(s: &SecretKey, q: usize, N: usize, p: usize) -> Self {
        let a = Polynomial::from_slice(&[1, 0, 0, 0]);
        let ring = ResidualQuotientRing::cyclotomic(N, q);

        let ek_1 = ring.add(
            &ring.mul(&a, &s.0).scale(-1),
            &ring.mul(&s.0, &s.0).scale(p as i64),
        );

        println!("EK: {}, {}", ek_1, a);
        Self(ek_1, a)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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

        assert_eq!(decoded.0, expected);
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

        assert_eq!(decoded.0, expected);
    }

    // WHAT WHY
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

        println!(
            "mul_cipher_and_decrypt: {}",
            session.mul_cipher_and_decrypt(c1.clone(), c2.clone()).0
        );

        println!(
            "mul_cipher_and_decrypt_decoded: {:?}",
            session
                .decode(session.mul_cipher_and_decrypt(c1.clone(), c2.clone()))
                .0
        );

        let c_mul = session.mul_cipher(c1, c2);

        let decrypted = session.decrypt(c_mul);
        println!("decrypted: {}", decrypted.0);
        println!(
            "huh: {}",
            decrypted
                .0
                .div_rem(Polynomial::from_slice(&[1, 0, 0, 0, 1]))
                .1
        );

        let expected = vec![Complex::new(9.0, 16.0), Complex::new(4.0, 1.0)];
        let pt = PlainText(expected);
        println!("expected: {}", session.encode(pt).0);

        let decoded = session.decode(decrypted);
        println!("decoded: {:?}", decoded.0);

        assert!(false);
    }
}
