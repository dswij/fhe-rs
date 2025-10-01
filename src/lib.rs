use poly::Polynomial;

extern crate nalgebra as nalg;

pub mod ckks;
pub mod dft;
pub mod poly;
mod session;

// const TEST_MOD: usize = 1024;
// const QUOTIENT_0: usize = 4096;
// const DELTA: usize = 16;
// const DELTA_L: usize = 3;
// const QUOTIENT_DELTA: usize = QUOTIENT_0 + DELTA.pow(DELTA_L as u32);
// const EK_CONSTANT: usize = QUOTIENT_DELTA;
// const N: usize = 1024;
//
// #[derive(Debug, Clone)]
// pub struct SecretKey(Polynomial<i64>);
//
// impl SecretKey {
//     fn generate_random() -> Self {
//         Self(gen_random_polynomial(1, 1))
//     }
// }
//
// #[derive(Debug, Clone)]
// pub struct PublicKey(Polynomial<i64>, Polynomial<i64>);
//
// impl PublicKey {
//     fn generate_from_secret_key(s: &SecretKey) -> Self {
//         let a = Polynomial::from_slice(&[0, 0, 0, 5]);
//         // let a = gen_random_polynomial(4, 4);
//         let e = gen_random_polynomial(1, 1);
//
//         let ring = ResidualQuotientRing::cyclotomic(N, TEST_MOD as i64);
//         let pk_1 = ring.add(&ring.mul(&a, &s.0), &e);
//
//         Self(pk_1.scale(-1), a)
//     }
// }
//
// #[derive(Debug, Clone)]
// struct CipherText(Polynomial<i64>, Polynomial<i64>);
//
// impl CipherText {
//     fn add(&self, rhs: &Self) -> Self {
//         let ring = ResidualQuotientRing::cyclotomic(N, TEST_MOD as i64);
//
//         Self(ring.add(&self.0, &rhs.0), ring.add(&self.1, &rhs.1))
//     }
//
//     fn mul(&self, rhs: &Self, ek: EvalKey) -> Self {
//         let c0 = self.0.clone() * rhs.0.clone();
//         let c1 = self.0.clone() * rhs.1.clone() + self.1.clone() * rhs.0.clone();
//
//         let scaler = TEST_MOD as i64;
//         let c2 = (self.1.clone() * rhs.1.clone()).scale(scaler as i64);
//
//         let ring = ResidualQuotientRing::cyclotomic(N, TEST_MOD as i64);
//         let c2 = (ring.mul(&ek.0, &c2), ring.mul(&ek.1, &c2));
//
//         Self(c0 + c2.0, c1 + c2.1)
//     }
// }
//
// #[derive(Debug, Clone)]
// struct Message(Polynomial<i64>);
//
// struct EvalKey(Polynomial<i64>, Polynomial<i64>);
//
// impl EvalKey {
//     fn generate_from_secret_key(s: &SecretKey) -> Self {
//         let a = gen_random_polynomial(1, 1);
//         let p = QUOTIENT_DELTA as i64;
//         let ring = ResidualQuotientRing {
//             ring: QuotientRing {
//                 quotient: Polynomial::from_slice(&[1, 0, 0, 0, 1]),
//             },
//             bound: TEST_MOD as i64,
//         };
//
//         let ek_1 = ring.add(
//             &ring.mul(&a, &s.0).scale(-1),
//             &ring.mul(&s.0, &s.0).scale(p),
//         );
//
//         Self(ek_1, a)
//     }
// }
//
// // TODO: figure out the bound for error and assert that
// pub fn encrypt(m: Message, pk: PublicKey) -> CipherText {
//     CipherText(m.0 + pk.0, pk.1)
// }
//
// pub fn decrypt(c: CipherText, sk: SecretKey) -> Message {
//     let ring = ResidualQuotientRing::cyclotomic(N, TEST_MOD as i64);
//
//     Message(c.0 + ring.mul(&c.1, &sk.0))
// }
//
// #[cfg(test)]
// mod tests {
//     use super::*;
//
//     use num::Complex;
//
//     #[test]
//     fn test_encrypt_decrypt() {
//         let sk = SecretKey::generate_random();
//         let pk = PublicKey::generate_from_secret_key(&sk);
//
//         let input = vec![Complex::new(3.0, 4.0), Complex::new(2.0, -1.0)];
//         let codec = ckks::codec::Codec::new(4, 64.0);
//         let msg = Message(codec.encode(input.clone()));
//
//         let cipher = encrypt(msg, pk);
//         let decrypted = decrypt(cipher, sk);
//
//         let decoded = codec.decode(decrypted.0);
//         dbg!(decoded);
//
//         assert!(false);
//     }
//
//     #[test]
//     fn test_addition() {
//         // let sk = SecretKey(Polynomial::from_slice(&[1, 0, 0, 100]));
//         let sk = SecretKey::generate_random();
//         let pk = PublicKey::generate_from_secret_key(&sk);
//
//         // cipher1
//         let input = vec![Complex::new(3.0, 4.0), Complex::new(2.0, -1.0)];
//         let codec = ckks::codec::Codec::new(4, 64.0);
//         let msg = Message(codec.encode(input.clone()));
//         let c1 = encrypt(msg, pk.clone());
//
//         // cipher2
//         let input = vec![Complex::new(3.0, 4.0), Complex::new(2.0, -1.0)];
//         let codec = ckks::codec::Codec::new(4, 64.0);
//         let msg = Message(codec.encode(input.clone()));
//         let c2 = encrypt(msg, pk);
//
//         let sum = c1.add(&c2);
//
//         let decrypted = decrypt(sum, sk);
//         let decoded = codec.decode(decrypted.0);
//
//         dbg!(decoded);
//
//         assert!(false);
//     }
//
//     #[test]
//     fn test_mult() {
//         let sk = SecretKey::generate_random();
//         let pk = PublicKey::generate_from_secret_key(&sk);
//
//         // cipher1
//         let input = vec![Complex::new(3.0, 4.0), Complex::new(2.0, -1.0)];
//         let codec = ckks::codec::Codec::new(4, 64.0);
//         let msg = Message(codec.encode(input.clone()));
//         let c1 = encrypt(msg, pk.clone());
//
//         // cipher2
//         let input = vec![Complex::new(2.0, 1.0), Complex::new(4.0, -1.0)];
//         let codec = ckks::codec::Codec::new(4, 64.0);
//         let msg = Message(codec.encode(input.clone()));
//         let c2 = encrypt(msg, pk);
//
//         let eval_key = EvalKey::generate_from_secret_key(&sk);
//         let multiplied = c1.mul(&c2, eval_key);
//
//         let decrypted = decrypt(multiplied, sk);
//         let decoded = codec.decode(decrypted.0);
//
//         dbg!(decoded);
//
//         assert!(false);
//     }
// }

// TODO: fix this
// pub struct PlainText(Polynomial<i64>);
// pub struct Message(Vec<u8>);
//
// struct Parameters {
//     N: isize,
//     q_0: isize,
//     L: isize,
// }
//
// pub const KEY_LENGTH: usize = 16;
// pub const BOUND: usize = 16; // the length bound for random small random polynomial
// pub const EK_CONSTANT: usize = 16;
// pub const Q: usize = 4096;
// pub const P: usize = 4096;
//
// pub fn encrypt(
//     m: Polynomial<i64>,
//     pk: (Polynomial<i64>, Polynomial<i64>),
// ) -> (Polynomial<i64>, Polynomial<i64>) {
//     (m + pk.0, pk.1)
// }
//
// pub fn decrypt(sk: Polynomial<i64>, cipher: (Polynomial<i64>, Polynomial<i64>)) -> Polynomial<i64> {
//     let mut plain_text = cipher.0 + cipher.1 * sk.clone();
//     // plain_text.coef_mod(Q.try_into().unwrap());
//     // plain_text.translate(-2048);
//
//     plain_text
// }
//
// /// Generates a random polynomial in R
// pub fn gen_secret_key(length: usize, bound: usize) -> Polynomial<i64> {
//     gen_random_polynomial(length, bound)
// }
//
// // Generates public key from a given secret key `sk`.
// pub fn gen_public_key(sk: Polynomial<i64>) -> (Polynomial<i64>, Polynomial<i64>) {
//     let e = gen_random_polynomial(KEY_LENGTH, 1);
//     let mut a = gen_random_polynomial(KEY_LENGTH, BOUND);
//
//     let mut pk = a.clone() * sk;
//     pk = pk.scale(-1);
//
//     // TODO: read more about the error `e` and it's bound.
//     pk = pk + e;
//     // pk.coef_mod(Q.try_into().unwrap());
//     // pk.translate(-2048);
//     // a.coef_mod(Q.try_into().unwrap());
//     // a.translate(-2048);
//
//     (pk, a)
// }
//
// fn gen_eval_key(secret_key: Polynomial<i64>) -> (Polynomial<i64>, Polynomial<i64>) {
//     let a_0 = gen_random_polynomial(KEY_LENGTH, BOUND);
//     let e_0 = gen_random_polynomial(KEY_LENGTH, BOUND);
//
//     let mut ek_0 = (a_0.clone() * secret_key.clone() + e_0).scale(-1)
//         + (secret_key.clone() * secret_key.clone()).scale(P as i64);
//
//     (ek_0, a_0)
// }
//
// fn gen_random_polynomial(length: usize, bound: usize) -> Polynomial<i64> {
//     let mut coefs = vec![0; length];
//
//     let bound = bound as i64;
//     for i in 0..length {
//         coefs[i] = rand::random_range(0..=bound);
//     }
//
//     Polynomial::from_slice(&coefs)
// }
//
// #[cfg(test)]
// mod tests {
//     use super::*;
//
//     use nalg::Complex;
//     use nalg::ComplexField;
//
//     #[test]
//     fn test_encrypt_decrypt() {
//         let sk = gen_secret_key(KEY_LENGTH, BOUND);
//         let pk = gen_public_key(sk.clone());
//
//         let input = vec![Complex::new(3.0, 4.0), Complex::new(2.0, -1.0)];
//         let codec = ckks::codec::Codec::new(4, 64.0);
//         let plain_text = codec.encode(input.clone());
//
//         let cipher = encrypt(plain_text.clone(), pk);
//
//         let decrypted = decrypt(sk, cipher);
//         let recovered: Vec<_> = codec
//             .decode(decrypted)
//             .into_iter()
//             .map(|e| e.round())
//             .collect();
//
//         assert_eq!(recovered, input);
//     }
//
//     #[test]
//     fn test_addition() {
//         let sk = gen_secret_key(KEY_LENGTH, BOUND);
//         let pk = gen_public_key(sk.clone());
//
//         let input_1 = vec![Complex::new(3.0, 4.0), Complex::new(2.0, 1.0)];
//         let input_2 = vec![Complex::new(6.0, 2.0), Complex::new(1.0, 0.0)];
//         let codec = ckks::codec::Codec::new(4, 64.0);
//         let pt_1 = codec.encode(input_1);
//         let pt_2 = codec.encode(input_2);
//         dbg!(&pt_1);
//         dbg!(&pt_2);
//
//         let cipher_1 = encrypt(pt_1, pk.clone());
//         let cipher_2 = encrypt(pt_2, pk);
//
//         let cipher = (cipher_1.0 + cipher_2.0, cipher_1.1 + cipher_2.1);
//         let decrypted = decrypt(sk, cipher);
//         let recovered: Vec<_> = codec
//             .decode(decrypted)
//             .into_iter()
//             .map(|e| e.round())
//             .collect();
//
//         let expected = vec![Complex::new(9.0, 6.0), Complex::new(3.0, 1.0)];
//         assert_eq!(recovered, expected);
//     }
//
//     #[test]
//     fn test_mult() {
//         let sk = gen_secret_key(KEY_LENGTH, BOUND);
//         let pk = gen_public_key(sk.clone());
//
//         let input_1 = vec![Complex::new(3.0, 0.0), Complex::new(2.0, 0.0)];
//         let input_2 = vec![Complex::new(6.0, 0.0), Complex::new(1.0, 0.0)];
//         let codec = ckks::codec::Codec::new(4, 64.0);
//         let pt_1 = codec.encode(input_1);
//         let pt_2 = codec.encode(input_2);
//         let eval_key = gen_eval_key(sk.clone());
//
//         let cipher_1 = encrypt(pt_1, pk.clone());
//         let cipher_2 = encrypt(pt_2, pk);
//
//         let cipher_mult = (
//             cipher_1.0.clone() * cipher_2.0.clone(),
//             cipher_1.0.clone() * cipher_2.1.clone() + cipher_1.1.clone() * cipher_2.0.clone(),
//         );
//
//         let cipher_mult_2 = cipher_1.1 * cipher_2.1;
//         let cipher_mult_relin = (
//             cipher_mult_2.clone() * eval_key.0,
//             cipher_mult_2.clone() * eval_key.1,
//         );
//
//         let decrypted = decrypt(sk.clone(), cipher_mult);
//         let recovered = codec.decode(decrypted).into_iter().map(|e| e.round());
//
//         let decrypted = decrypt(sk.clone(), cipher_mult_relin);
//         let recovered2: Vec<_> = codec
//             .decode(decrypted)
//             .into_iter()
//             .zip(recovered)
//             .map(|(e1, e2)| (e1.round().re + e2.re) % 256.0)
//             .collect();
//         dbg!(&recovered2);
//         assert!(false);
//         // let expected = vec![Complex::new(18.0, 0.0), Complex::new(2.0, 0.0)];
//         // assert_eq!(recovered2, expected);
//     }
// }
