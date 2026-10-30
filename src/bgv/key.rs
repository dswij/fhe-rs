use crate::poly::gen_random_polynomial;
use crate::poly::Polynomial;
use crate::sample::uniform_R2;

use crate::bgv::{BGVParams, KeySwitch};

#[derive(Clone, PartialEq)]
pub struct SecretKey(pub Polynomial<i64>);

impl SecretKey {
    pub fn new(p: Polynomial<i64>) -> Self {
        Self(p)
    }

    pub fn generate_random(len: usize) -> Self {
        Self::new(uniform_R2(len))
    }
}

#[derive(Debug, Clone)]
pub struct PublicKey(pub Polynomial<i64>, pub Polynomial<i64>);

impl PublicKey {
    // PK1 = [−1(a · SK + t · e)]_qL
    // PK2 = a
    pub fn generate_from_secret_key(s: &SecretKey, params: &BGVParams) -> Self {
        let ring = params.ring_at(params.top_level());
        let a = gen_random_polynomial(params.n, params.top_modulus());
        let e = gen_random_polynomial(params.n, 2);

        let mut pk_1 = ring
            .add(&ring.mul(&a, &s.0), &e.scale(params.t as i64))
            .scale(-1);
        pk_1.coef_mod(params.top_modulus() as i64);

        Self(pk_1, a)
    }
}

// Key-switching key from `s^2` to `s`, generated at the top level. Its shape depends on
// [`KeySwitch`].
#[derive(Debug, Clone)]
pub enum EvalKey {
    // One pair per base-`w` digit: `( -(a_i * s + t * e_i) + w^i * s^2 , a_i )`, all mod `Q_L`,
    // so `evk_i.0 + evk_i.1 * s = w^i * s^2 - t * e_i`.
    Digits(Vec<(Polynomial<i64>, Polynomial<i64>)>),

    // A single pair `( -(a * s + t * e) + P * s^2 , a )` modulo `P * Q_L`.
    ModulusExtension(Polynomial<i64>, Polynomial<i64>),
}

impl EvalKey {
    /// The error is a multiple of `t`, as in the public key, so it does not disturb the plaintext.
    pub fn generate_from_secret_key(s: &SecretKey, params: &BGVParams) -> Self {
        match params.key_switch {
            KeySwitch::Digits { w } => Self::generate_digits(s, params, w),
            KeySwitch::ModulusExtension { p } => Self::generate_special(s, params, p),
        }
    }

    /// The key for ciphertexts at `level`. `Q_level` divides `Q_L`, so reducing every coefficient
    /// keeps `k0 + k1 * s` equal to the same value modulo the smaller modulus.
    pub fn at_level(&self, params: &BGVParams, level: usize) -> Self {
        let reduce = |p: &Polynomial<i64>, m: usize| {
            let mut p = p.clone();
            p.coef_mod(m as i64);
            p
        };

        match (self, params.key_switch) {
            (Self::Digits(pairs), KeySwitch::Digits { w }) => {
                let q = params.modulus(level);
                Self::Digits(
                    pairs
                        .iter()
                        .take(params.num_digits(level, w))
                        .map(|(k0, k1)| (reduce(k0, q), reduce(k1, q)))
                        .collect(),
                )
            }
            (Self::ModulusExtension(k0, k1), KeySwitch::ModulusExtension { p }) => {
                let pq = p * params.modulus(level);
                Self::ModulusExtension(reduce(k0, pq), reduce(k1, pq))
            }
            (_, switch) => panic!("evaluation key does not match the key switch {switch:?}"),
        }
    }

    fn generate_digits(s: &SecretKey, params: &BGVParams, w: usize) -> Self {
        let top = params.top_level();
        let q = params.top_modulus();
        let ring = params.ring_at(top);
        let s2 = ring.mul(&s.0, &s.0);

        let pairs = (0..params.num_digits(top, w))
            .map(|i| {
                let a = gen_random_polynomial(params.n, q);
                let e = gen_random_polynomial(params.n, 2).scale(params.t as i64);

                let masked = ring.add(&ring.mul(&a, &s.0), &e).scale(-1);
                let wi = (w as i64).pow(i as u32);
                (ring.add(&masked, &s2.clone().scale(wi)), a)
            })
            .collect();

        Self::Digits(pairs)
    }

    fn generate_special(s: &SecretKey, params: &BGVParams, p: usize) -> Self {
        assert_eq!(
            p % params.t,
            1,
            "special modulus needs P = 1 (mod t), got P = {p}, t = {}",
            params.t
        );

        let ring = params.ring_pq_at(params.top_level(), p);
        let a = gen_random_polynomial(params.n, p * params.top_modulus());
        let e = gen_random_polynomial(params.n, 2).scale(params.t as i64);

        let masked = ring.add(&ring.mul(&a, &s.0), &e).scale(-1);
        let k0 = ring.add(&masked, &ring.mul(&s.0, &s.0).scale(p as i64));

        Self::ModulusExtension(k0, a)
    }
}
