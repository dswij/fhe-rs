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
