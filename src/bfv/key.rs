use crate::poly::gen_random_polynomial;
use crate::poly::ring::*;
use crate::poly::Polynomial;

use crate::bfv::BFVParams;

#[derive(Clone, PartialEq)]
pub struct SecretKey(pub Polynomial<i64>);

impl SecretKey {
    pub fn new(p: Polynomial<i64>) -> Self {
        Self(p)
    }
}

#[derive(Debug, Clone)]
pub struct PublicKey(pub Polynomial<i64>, pub Polynomial<i64>);

impl PublicKey {
    pub fn generate_from_secret_key(s: &SecretKey, params: BFVParams) -> Self {
        let a = gen_random_polynomial(params.N, params.Q);
        let e = gen_random_polynomial(params.N, 2);

        let ring = ResidualQuotientRing::cyclotomic(params.N, params.Q);
        let pk_1 = ring.add(&ring.mul(&a, &s.0), &e);

        Self(pk_1.scale(-1), a)
    }
}

#[derive(Debug, Clone)]
pub struct EvalKey(pub Polynomial<i64>, pub Polynomial<i64>);

impl EvalKey {
    pub fn generate_from_secret_key(s: &SecretKey, params: BFVParams) -> Self {
        let a = gen_random_polynomial(params.N, params.P * params.Q);
        let e = gen_random_polynomial(params.N, 2);
        let ring = ResidualQuotientRing::cyclotomic(params.N, params.P * params.Q);

        let masked = ring.add(&ring.mul(&a, &s.0), &e).scale(-1);
        let ek = ring.add(&masked, &ring.mul(&s.0, &s.0).scale(params.P as i64));

        Self(ek, a)
    }
}
