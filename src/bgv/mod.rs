//! BGV Scheme
//!
//! Implementation reference: https://faculty.kfupm.edu.sa/coe/mfelemban/SEC595/References/Introduction%20to%20the%20BGV%20FHE%20Scheme.pdf
pub mod chain;
pub mod key;

use crate::poly::gen_random_polynomial;
use crate::poly::ring::*;
use crate::poly::Polynomial;
use crate::sample::uniform_R2;

use chain::ModulusChain;
use key::*;

// TODO: write a sample app for this?
// TODO: figure out the noise calculation for this, especially EvalMul

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum KeySwitch {
    // Original BGV, relin by digit decomposition.
    Digits { w: usize },
    // GHS Variant, for Mul, using the modulo P * Q when relin instead of modulo Q.
    ModulusExtension { p: usize },
}

#[derive(Debug, Clone)]
pub struct BGVParams {
    pub n: usize,
    pub t: usize,
    pub chain: ModulusChain,
    pub key_switch: KeySwitch,
}

impl BGVParams {
    pub fn top_level(&self) -> usize {
        self.chain.top_level()
    }

    pub fn modulus(&self, level: usize) -> usize {
        self.chain.modulus(level)
    }

    pub fn top_modulus(&self) -> usize {
        self.chain.top()
    }

    fn ring_at(&self, level: usize) -> ResidualQuotientRing {
        ResidualQuotientRing::cyclotomic(self.n, self.modulus(level))
    }

    fn ring_pq_at(&self, level: usize, p: usize) -> ResidualQuotientRing {
        ResidualQuotientRing::cyclotomic(self.n, p * self.modulus(level))
    }

    // Number of base-`w` digits needed to represent a value in `[0, Q_level)`.
    pub fn num_digits(&self, level: usize, w: usize) -> usize {
        let (mut digits, mut x) = (0, 1usize);
        while x < self.modulus(level) {
            x *= w;
            digits += 1;
        }
        digits
    }
}

// Implement Eval()
pub struct Evaluator {
    params: BGVParams,
}

impl Evaluator {
    // EvalAdd
    pub fn add(&self, lhs: CipherText, rhs: CipherText) -> CipherText {
        assert_eq!(lhs.2, rhs.2, "ciphertexts are at different levels");
        let ring = self.params.ring_at(lhs.2);

        CipherText(ring.add(&lhs.0, &rhs.0), ring.add(&lhs.1, &rhs.1), lhs.2)
    }

    // EvalMul.
    pub fn mul(&self, lhs: CipherText, rhs: CipherText, ek: EvalKey) -> CipherText {
        assert_eq!(lhs.2, rhs.2, "ciphertexts are at different levels");
        let level = lhs.2;
        let ring = self.params.ring_at(level);
        let c0 = ring.mul(&lhs.0, &rhs.0);
        let c1 = ring.add(&ring.mul(&lhs.0, &rhs.1), &ring.mul(&lhs.1, &rhs.0));
        let c2 = ring.mul(&lhs.1, &rhs.1);

        match (self.params.key_switch, ek.at_level(&self.params, level)) {
            (KeySwitch::Digits { w }, EvalKey::Digits(pairs)) => {
                self.relinearize_digits(&c0, &c1, &c2, &pairs, w as i64, level)
            }
            (KeySwitch::ModulusExtension { p }, EvalKey::ModulusExtension(k0, k1)) => {
                self.relinearize_special(&c0, &c1, &c2, &k0, &k1, p as i64, level)
            }
            (switch, _) => {
                panic!("unexpected switch or switch/evalkey combination. switch: {switch:?}");
            }
        }
    }

    // BGV key switching: write c2 = sum_i d_i * w^i with small digits d_i, then add
    // sum_i d_i * evk_i. No division, so there is nothing to correct mod t.
    fn relinearize_digits(
        &self,
        c0: &Polynomial<i64>,
        c1: &Polynomial<i64>,
        c2: &Polynomial<i64>,
        pairs: &[(Polynomial<i64>, Polynomial<i64>)],
        w: i64,
        level: usize,
    ) -> CipherText {
        let ring = self.params.ring_at(level);
        let (mut c0, mut c1) = (c0.clone(), c1.clone());
        for (i, (k0, k1)) in pairs.iter().enumerate() {
            let wi = w.pow(i as u32);
            let digit = Polynomial::new(c2.coef.iter().map(|x| (x / wi) % w).collect());

            c0 = ring.add(&c0, &ring.mul(&digit, k0));
            c1 = ring.add(&c1, &ring.mul(&digit, k1));
        }

        CipherText(c0, c1, level)
    }

    // GHS key switching: the key is `(-(a*s + t*e) + P*s^2, a)` mod `P*Q`. Compute
    // `(P*c0 + k0*c2, P*c1 + k1*c2)` mod `P*Q`, then scale back down to `Q`. The added noise
    // `t * c2 * e` is divided by `P` by the scaling.
    fn relinearize_special(
        &self,
        c0: &Polynomial<i64>,
        c1: &Polynomial<i64>,
        c2: &Polynomial<i64>,
        k0: &Polynomial<i64>,
        k1: &Polynomial<i64>,
        p: i64,
        level: usize,
    ) -> CipherText {
        let ring_pq = self.params.ring_pq_at(level, p as usize);
        let x0 = ring_pq.add(&c0.clone().scale(p), &ring_pq.mul(k0, c2));
        let x1 = ring_pq.add(&c1.clone().scale(p), &ring_pq.mul(k1, c2));

        let (pq, q, t) = (
            p * self.params.modulus(level) as i64,
            self.params.modulus(level) as i64,
            self.params.t as i64,
        );
        CipherText(scale(&x0, pq, q, t), scale(&x1, pq, q, t), level)
    }

    // Q_i -> Q_{i-1}: divide by the dropped prime `q_i`. The noise goes
    // from `v` to about `v / q_i + (t / 2) * (1 + ||s||_1)`.
    pub fn mod_switch(&self, cipher: &CipherText) -> CipherText {
        let level = cipher.level();
        assert!(level > 0, "already at the lowest level");

        let t = self.params.t;
        let dropped = self.params.chain.prime(level);
        // `Q_level = Q_{level-1} * dropped`, so the two moduli agree mod `t` exactly when
        // `dropped = 1 (mod t)`
        //
        // This assertion is what keeps the plaintext intact(Lemma 4 of the BGV paper).
        // TODO: maybe we don't need to assert in release, only for dev profile.
        assert_eq!(
            dropped % t,
            1,
            "cannot switch from level {level}: q_{level} = {dropped} must be 1 (mod t = {t})"
        );

        let (q, q2) = (
            self.params.modulus(level) as i64,
            self.params.modulus(level - 1) as i64,
        );
        CipherText(
            scale(&cipher.0, q, q2, t as i64),
            scale(&cipher.1, q, q2, t as i64),
            level - 1,
        )
    }
}

// `Scale(x, q, q2, t)`: per coefficient, the integer closest to `v * q2 / q` that is congruent
// to `v` mod `t`, reduced mod `q2`. `v` is in `[0, q)`.
//
// Here's roughly what is done:
//
// Let y be any integer that is congruent to v. i.e. y mod t = v mod t = r.
// For such y, we can write y = r + t * j.
//
// We want to find an integer that is closest to the value of v / q, but maintaining the
// congruence. That is, minimizing the distance between v / q and y. Replace v with v * q2.
//
// `|v*q2/q - r - t*j| = t * |(v*q2 - q*r) / (q*t) - j|`, so choose the smallest j.
//
// Also note that `round(a / b) = floor((2a + b) / (2b))`, which is always the trick to do
// rounding... without round function.
//
// The products reach `q^2`, so work in i128.
fn scale(x: &Polynomial<i64>, q: i64, q2: i64, t: i64) -> Polynomial<i64> {
    let (q, q2, t) = (q as i128, q2 as i128, t as i128);
    Polynomial::new(
        x.coef
            .iter()
            .map(|v| {
                let v = *v as i128;
                let r = v.rem_euclid(t);
                let (num, den) = (v * q2 - q * r, q * t);
                let j = (2 * num + den).div_euclid(2 * den);
                (r + t * j).rem_euclid(q2) as i64
            })
            .collect(),
    )
}

pub struct Encrypt {
    params: BGVParams,
    pk: PublicKey,
}

impl Encrypt {
    pub fn encrypt(&self, m: Message) -> CipherText {
        let level = self.params.top_level();
        let ring = self.params.ring_at(level);
        let n = self.params.n;
        let t = self.params.t as i64;
        let mut u = uniform_R2(n);
        while u.is_zero() || u.coef.iter().all(|e| *e == 0) {
            u = uniform_R2(n);
        }

        let e1 = gen_random_polynomial(n, 2).scale(t);
        let e2 = gen_random_polynomial(n, 2).scale(t);
        let c0 = ring.add(&ring.add(&ring.mul(&self.pk.0, &u), &e1), &m.0.clone());
        let c1 = ring.add(&ring.mul(&self.pk.1, &u), &e2);

        CipherText(c0, c1, level)
    }
}

pub struct Decrypt {
    params: BGVParams,
    sk: SecretKey,
}

impl Decrypt {
    pub fn decrypt(&self, c: &CipherText) -> Message {
        let ring = self.params.ring_at(c.2);
        let (q, t) = (self.params.modulus(c.2) as i64, self.params.t as i64);

        let mut v = ring.add(&c.0, &ring.mul(&c.1, &self.sk.0)); // reduced mod Q_level
        v.coef.iter_mut().for_each(|x| {
            if *x > q / 2 {
                *x -= q
            } // center lift
            *x = x.rem_euclid(t);
        });
        Message(v)
    }
}

/// `(c0, c1, level)`: the polynomials are modulo `Q_level`.
#[derive(Debug, Clone)]
pub struct CipherText(Polynomial<i64>, Polynomial<i64>, usize);

impl CipherText {
    pub fn level(&self) -> usize {
        self.2
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Message(Polynomial<i64>);

#[cfg(test)]
mod tests {
    use super::*;
    use tests::utils::*;

    const TRIAL_COUNT: usize = 30;

    #[test]
    fn test_encrypt_decrypt() {
        let scheme = test_scheme();

        let input = Polynomial::from_slice(&[5, 4]);
        let msg = Message(input);
        let encrypted = scheme.enc.encrypt(msg);
        let decrypted = scheme.dec.decrypt(&encrypted);

        assert_eq!(decrypted.0, Polynomial::from_slice(&[5, 4, 0, 0]));
    }

    // A ciphertext must be masked by `pk * u`: its coefficients should be spread over Z_Q,
    // not just `m` plus a few units of noise.
    #[test]
    fn test_ciphertext_is_masked() {
        let q = test_params().top_modulus() as i64;
        let msg = Polynomial::from_slice(&[5, 4, 0, 0]);

        let (mut c0_masked, mut c1_masked) = (0, 0);
        for _ in 0..50 {
            let scheme = test_scheme();
            let c = scheme.enc.encrypt(Message(msg.clone()));

            let mask = c.0.clone() - msg.clone();
            c0_masked += (max_centered(&mask, q) > q / 16) as usize;
            c1_masked += (max_centered(&c.1, q) > q / 16) as usize;
        }

        // The mask vanishes when the ternary `u` is zero (about 1 in 81), so allow a few misses.
        assert!(c0_masked >= 40, "c0 was masked in only {c0_masked} of 50");
        assert!(c1_masked >= 40, "c1 was masked in only {c1_masked} of 50");
    }

    // Decrypting under a different secret key must not reveal the message.
    #[test]
    fn test_wrong_secret_key_cannot_decrypt() {
        // A wrong key decrypts to an effectively random polynomial, so it matches `msg` with
        // probability 1/t^N. Make this longer so we virtually won't run into this case.
        let params = BGVParams {
            n: 16,
            ..test_params()
        };
        let msg = Polynomial::from_slice(&[5, 4, 3, 2, 1, 7, 6, 5, 4, 3, 2, 1, 7, 6, 5, 4]);

        for _ in 0..TRIAL_COUNT {
            let (alice, mut bob) = (
                scheme_with_params(params.clone()),
                scheme_with_params(params.clone()),
            );
            // Ensure that we'll never hit the case where it's the same SK generated.
            while bob.dec.sk == alice.dec.sk {
                bob = scheme_with_params(params.clone());
            }
            let c = alice.enc.encrypt(Message(msg.clone()));
            assert!(bob.dec.decrypt(&c).0 != msg);
        }
    }

    #[test]
    fn test_add() {
        // (t, Q, lhs, rhs, expected)
        let cases: Vec<(usize, usize, Vec<i64>, Vec<i64>, Vec<i64>)> = vec![
            (8, 525641, vec![1, 2], vec![1, 2], vec![2, 4, 0, 0]),
            // wraparound: (7 + 3) % 2 = 2
            (8, 525641, vec![7, 0], vec![3, 0], vec![2, 0, 0, 0]),
        ];

        for (t, q, a, b, expected) in cases {
            let scheme = scheme_with(t, q);
            let expected = Polynomial::from_slice(&expected);
            let ma = Message(Polynomial::from_slice(&a));
            let mb = Message(Polynomial::from_slice(&b));

            for _ in 0..TRIAL_COUNT {
                let ca = scheme.enc.encrypt(ma.clone());
                let cb = scheme.enc.encrypt(mb.clone());

                let sum = scheme.dec.decrypt(&scheme.eval.add(ca, cb));
                assert_eq!(sum.0, expected, "t = {t}, Q = {q}, {a:?} * {b:?}");
            }
        }
    }

    #[test]
    fn test_mul() {
        // (t, Q, lhs, rhs, expected)
        let cases: Vec<(usize, usize, Vec<i64>, Vec<i64>, Vec<i64>)> = vec![
            (8, 525641, vec![1, 2], vec![1, 2], vec![1, 4, 4, 0]),
            (8, 525641, vec![1, 3], vec![1, 3], vec![1, 6, 1, 0]),
            // 3 * 3 = 9, wraps to 1 mod t
            (8, 525641, vec![3], vec![3], vec![1, 0, 0, 0]),
            // x^3 * x = x^4 = -1 in the negacyclic ring, i.e. t - 1
            (8, 525641, vec![0, 0, 0, 1], vec![0, 1], vec![7, 0, 0, 0]),
            (17, 525641, vec![1, 2], vec![3, 1], vec![3, 7, 2, 0]),
            (
                17,
                525641,
                vec![5, 9, 2, 16],
                vec![7, 0, 4, 1],
                vec![1, 14, 1, 0],
            ),
        ];

        for key_switch in key_switches() {
            for (t, q, a, b, expected) in &cases {
                let (t, q) = (*t, *q);
                let scheme = scheme_with_key_switch(t, q, key_switch);
                let expected = Polynomial::from_slice(expected);
                let ma = Message(Polynomial::from_slice(a));
                let mb = Message(Polynomial::from_slice(b));

                for _ in 0..TRIAL_COUNT {
                    let ca = scheme.enc.encrypt(ma.clone());
                    let cb = scheme.enc.encrypt(mb.clone());

                    let product = scheme
                        .dec
                        .decrypt(&scheme.eval.mul(ca, cb, scheme.ek.clone()));

                    assert_eq!(
                        product.0, expected,
                        "{key_switch:?}, t = {t}, Q = {q}, {a:?} * {b:?}"
                    );
                }
            }
        }
    }

    // The special-modulus variant only preserves the plaintext when `P = 1 (mod t)`.
    #[test]
    #[should_panic(expected = "P = 1 (mod t)")]
    fn test_special_modulus_requires_p_congruent_one() {
        // 1025 = 1 (mod 8) but 1025 = 5 (mod 17).
        scheme_with_key_switch(17, 525641, KeySwitch::ModulusExtension { p: 1025 });
    }

    #[test]
    fn test_random_messages() {
        let (t, q) = (17usize, 525641usize);
        for key_switch in key_switches() {
            for _ in 0..TRIAL_COUNT {
                let scheme = scheme_with_key_switch(t, q, key_switch);
                let rand_msg =
                    || -> Vec<i64> { (0..4).map(|_| rand::random_range(0..t as i64)).collect() };
                let (a, b) = (rand_msg(), rand_msg());

                let ca = scheme.enc.encrypt(Message(Polynomial::from_slice(&a)));
                let cb = scheme.enc.encrypt(Message(Polynomial::from_slice(&b)));

                let sum = scheme.dec.decrypt(&scheme.eval.add(ca.clone(), cb.clone()));
                let expected_sum: Vec<i64> =
                    a.iter().zip(&b).map(|(x, y)| (x + y) % t as i64).collect();
                assert_eq!(
                    sum.0,
                    Polynomial::from_slice(&expected_sum),
                    "{key_switch:?}: {a:?} + {b:?}"
                );

                let product = scheme.dec.decrypt(&scheme.eval.mul(ca, cb, scheme.ek));
                assert_eq!(
                    product.0,
                    plain_mul(&a, &b, 4, t as i64),
                    "{key_switch:?}: {a:?} * {b:?}"
                );
            }
        }
    }

    // `mod_switch` moves the ciphertext down one level without changing what it decrypts to.
    #[test]
    fn test_mod_switch_preserves_message() {
        for t in [8usize, 17] {
            let chain = ModulusChain::generate(t, 4, &[10, 9, 9]);
            let scheme = scheme_with_chain(t, chain, test_params().key_switch);

            for _ in 0..TRIAL_COUNT {
                let msg = random_message(t);
                let c = scheme.enc.encrypt(Message(Polynomial::from_slice(&msg)));
                assert_eq!(c.level(), 2);

                let expected = Polynomial::from_slice(&msg);
                let c = scheme.eval.mod_switch(&c);
                assert_eq!(c.level(), 1);
                assert_eq!(scheme.dec.decrypt(&c).0, expected, "t = {t}, level 1");

                let c = scheme.eval.mod_switch(&c);
                assert_eq!(c.level(), 0);
                assert_eq!(scheme.dec.decrypt(&c).0, expected, "t = {t}, level 0");
            }
        }
    }

    // After a switch the noise is at most `v * q' / q + (t / 2) * (1 + ||s||_1)`, with
    // `||s||_1 <= N`. Check the bound on a product, whose noise is far above that floor.
    #[test]
    fn test_mod_switch_noise_bound() {
        let (t, n) = (8usize, 4usize);
        let chain = ModulusChain::generate(t, n, &[10, 9, 9]);
        let dropped = chain.prime(2) as f64;
        let scheme = scheme_with_chain(t, chain, test_params().key_switch);
        let floor = (t as f64 / 2.0) * (1.0 + n as f64);

        for _ in 0..TRIAL_COUNT {
            let ca = scheme
                .enc
                .encrypt(Message(Polynomial::from_slice(&random_message(t))));
            let cb = scheme
                .enc
                .encrypt(Message(Polynomial::from_slice(&random_message(t))));
            let product = scheme.eval.mul(ca, cb, scheme.ek.clone());

            let before = phase_norm(&scheme, &product) as f64;
            let after = phase_norm(&scheme, &scheme.eval.mod_switch(&product)) as f64;

            // +1 for the integer rounding of the centered norm.
            assert!(
                after <= before / dropped + floor + 1.0,
                "noise {before} became {after}, expected at most {}",
                before / dropped + floor + 1.0
            );
            assert!(
                after < before / 10.0,
                "noise barely shrank: {before} -> {after}"
            );
        }
    }

    // Two multiplications in a row need a switch between them: (a * b) * c at three levels.
    #[test]
    fn test_depth_two_with_mod_switch() {
        let t = 17usize;
        let chain = ModulusChain::generate(t, 4, &[10, 9, 9]);
        let scheme = scheme_with_chain(t, chain, KeySwitch::Digits { w: 128 });

        for _ in 0..TRIAL_COUNT {
            let (a, b, c) = (random_message(t), random_message(t), random_message(t));
            let enc = |m: &[i64]| scheme.enc.encrypt(Message(Polynomial::from_slice(m)));

            let ab = scheme.eval.mul(enc(&a), enc(&b), scheme.ek.clone());
            let ab = scheme.eval.mod_switch(&ab);
            // `c` has to come down to the same level before it can be multiplied.
            let c_ct = scheme.eval.mod_switch(&enc(&c));
            let abc = scheme.eval.mul(ab, c_ct, scheme.ek.clone());
            let abc = scheme.eval.mod_switch(&abc);

            assert_eq!(abc.level(), 0);
            let expected = plain_mul(&plain_mul(&a, &b, 4, t as i64).coef, &c, 4, t as i64);
            assert_eq!(
                scheme.dec.decrypt(&abc).0,
                expected,
                "{a:?} * {b:?} * {c:?}"
            );
        }
    }

    // The same with the special-modulus key switch. `P * Q_1` must stay inside i64 products, so
    // the chain has two levels.
    #[test]
    fn test_mod_switch_after_special_modulus_mul() {
        let t = 17usize;
        let chain = ModulusChain::generate(t, 4, &[10, 9]);
        let scheme = scheme_with_chain(t, chain, KeySwitch::ModulusExtension { p: 1089 });

        for _ in 0..TRIAL_COUNT {
            let (a, b) = (random_message(t), random_message(t));
            let enc = |m: &[i64]| scheme.enc.encrypt(Message(Polynomial::from_slice(m)));

            let ab = scheme.eval.mul(enc(&a), enc(&b), scheme.ek.clone());
            let ab = scheme.eval.mod_switch(&ab);

            assert_eq!(
                scheme.dec.decrypt(&ab).0,
                plain_mul(&a, &b, 4, t as i64),
                "{a:?} * {b:?}"
            );
        }
    }

    #[test]
    #[should_panic(expected = "already at the lowest level")]
    fn test_mod_switch_below_level_zero_panics() {
        let chain = ModulusChain::generate(8, 4, &[10, 9]);
        let scheme = scheme_with_chain(8, chain, test_params().key_switch);
        let c = scheme.enc.encrypt(Message(Polynomial::from_slice(&[1])));

        scheme.eval.mod_switch(&scheme.eval.mod_switch(&c));
    }

    // A dropped prime that is not 1 mod t would change the plaintext, so it is refused.
    #[test]
    #[should_panic(expected = "must be 1 (mod t")]
    fn test_mod_switch_requires_prime_congruent_one() {
        // 1009 = 1 (mod 8) but 461 = 5 (mod 8).
        let chain = ModulusChain::new(vec![1009, 461]);
        let scheme = scheme_with_chain(8, chain, test_params().key_switch);
        let c = scheme.enc.encrypt(Message(Polynomial::from_slice(&[1])));

        scheme.eval.mod_switch(&c);
    }

    pub mod utils {
        use super::*;

        pub struct Scheme {
            pub enc: Encrypt,
            pub dec: Decrypt,
            pub eval: Evaluator,
            pub ek: EvalKey,
        }

        pub fn test_params() -> BGVParams {
            BGVParams {
                n: 4,
                t: 8,
                chain: ModulusChain::new(vec![525641]),
                key_switch: KeySwitch::Digits { w: 128 },
            }
        }

        // Both key-switching modes. `P = 1 (mod 8)` and `P = 1 (mod 17)`, so one `P` serves
        // every plaintext modulus used in the tests.
        pub fn key_switches() -> [KeySwitch; 2] {
            [
                KeySwitch::Digits { w: 128 },
                KeySwitch::ModulusExtension { p: 1089 },
            ]
        }

        pub fn scheme_with(t: usize, q: usize) -> Scheme {
            scheme_with_key_switch(t, q, test_params().key_switch)
        }

        pub fn scheme_with_key_switch(t: usize, q: usize, key_switch: KeySwitch) -> Scheme {
            scheme_with_chain(t, ModulusChain::new(vec![q]), key_switch)
        }

        pub fn scheme_with_chain(t: usize, chain: ModulusChain, key_switch: KeySwitch) -> Scheme {
            scheme_with_params(BGVParams {
                t,
                chain,
                key_switch,
                ..test_params()
            })
        }

        pub fn scheme_with_params(params: BGVParams) -> Scheme {
            let sk = SecretKey::generate_random(params.n);
            let pk = PublicKey::generate_from_secret_key(&sk, &params);
            let ek = EvalKey::generate_from_secret_key(&sk, &params);

            Scheme {
                enc: Encrypt {
                    params: params.clone(),
                    pk,
                },
                dec: Decrypt {
                    params: params.clone(),
                    sk,
                },
                eval: Evaluator { params },
                ek,
            }
        }

        pub fn test_scheme() -> Scheme {
            let params = test_params();
            scheme_with(params.t, params.top_modulus())
        }

        // Largest coefficient of the centered phase `c0 + c1 * s` of `c`: the message plus noise.
        pub fn phase_norm(scheme: &Scheme, c: &CipherText) -> i64 {
            let ring = scheme.dec.params.ring_at(c.2);
            let phase = ring.add(&c.0, &ring.mul(&c.1, &scheme.dec.sk.0));

            max_centered(&phase, scheme.dec.params.modulus(c.2) as i64)
        }

        pub fn random_message(t: usize) -> Vec<i64> {
            (0..4).map(|_| rand::random_range(0..t as i64)).collect()
        }

        // Reference product in Z_t[X]/(X^n + 1).
        pub fn plain_mul(a: &[i64], b: &[i64], n: usize, t: i64) -> Polynomial<i64> {
            let mut c = vec![0i64; n];
            for (i, x) in a.iter().enumerate() {
                for (j, y) in b.iter().enumerate() {
                    let k = i + j;
                    if k < n {
                        c[k] += x * y;
                    } else {
                        c[k - n] -= x * y;
                    }
                }
            }
            c.iter_mut().for_each(|e| *e = e.rem_euclid(t));
            Polynomial::new(c)
        }

        // Largest coefficient magnitude of `p` after centering mod `q`.
        pub fn max_centered(p: &Polynomial<i64>, q: i64) -> i64 {
            p.coef
                .iter()
                .map(|x| {
                    let x = x.rem_euclid(q);
                    x.min(q - x)
                })
                .max()
                .unwrap_or(0)
        }
    }
}
