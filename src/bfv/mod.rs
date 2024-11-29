//! BFV Scheme
//!
//! Implementation reference: https://faculty.kfupm.edu.sa/coe/mfelemban/SEC595/References/Introduction%20to%20the%20BFV%20FHE%20Scheme.pdf

pub mod key;

use crate::poly::gen_random_polynomial;
use crate::poly::ring::*;
use crate::poly::Polynomial;
use crate::sample::uniform_R2;

use key::*;

#[derive(Debug, Clone, Copy)]
pub struct BFVParams {
    N: usize,
    t: usize,
    Q: usize,
    P: usize,
}

impl BFVParams {
    fn delta(&self) -> i64 {
        (self.Q / self.t) as i64
    }

    fn ring_q(&self) -> ResidualQuotientRing {
        ResidualQuotientRing::cyclotomic(self.N, self.Q)
    }

    fn ring_pq(&self) -> ResidualQuotientRing {
        ResidualQuotientRing::cyclotomic(self.N, self.P * self.Q)
    }
}

// Implement Eval()
pub struct Evaluator {
    params: BFVParams,
}

impl Evaluator {
    // EvalAdd
    pub fn add(&self, lhs: CipherText, rhs: CipherText) -> CipherText {
        let ring = self.params.ring_q();

        CipherText(ring.add(&lhs.0, &rhs.0), ring.add(&lhs.1, &rhs.1))
    }

    // EvalMul
    pub fn mul(&self, lhs: CipherText, rhs: CipherText, ek: EvalKey) -> CipherText {
        let (t, q) = (self.params.t as i64, self.params.Q as i64);

        let c0 = &lhs.0 * &rhs.0;
        let c0 = Polynomial::new(
            c0.coef
                .iter()
                .map(|m| (2 * t * m + q).div_euclid(2 * q).rem_euclid(q))
                .collect(),
        );
        let c1 = &lhs.0 * &rhs.1 + &lhs.1 * &rhs.0;
        let c1 = Polynomial::new(
            c1.coef
                .iter()
                .map(|m| (2 * t * m + q).div_euclid(2 * q).rem_euclid(q))
                .collect(),
        );

        let c2 = &lhs.1 * &rhs.1;
        let c2 = Polynomial::new(
            c2.coef
                .iter()
                .map(|m| (2 * t * m + q).div_euclid(2 * q).rem_euclid(q))
                .collect(),
        );

        let ring_pq = self.params.ring_pq();
        let (p, q) = (self.params.P as i64, self.params.Q as i64);
        let mod_down = |x: Polynomial<i64>| {
            Polynomial::new(
                x.coef
                    .iter()
                    .map(|v| ((v + p / 2).div_euclid(p)).rem_euclid(q))
                    .collect(),
            )
        };

        let r0 = mod_down(ring_pq.mul(&ek.0, &c2));
        let r1 = mod_down(ring_pq.mul(&ek.1, &c2));

        let ring = self.params.ring_q();
        let c0 = ring.add(&c0, &r0);
        let c1 = ring.add(&c1, &r1);

        CipherText(c0, c1)
    }
}

pub struct Encrypt {
    params: BFVParams,
    pk: PublicKey,
}

impl Encrypt {
    pub fn encrypt(&self, m: Message) -> CipherText {
        let ring = self.params.ring_q();
        let n = self.params.N;
        let u = uniform_R2(n);
        let e1 = gen_random_polynomial(n, 2);
        let e2 = gen_random_polynomial(n, 2);
        let delta = self.params.delta();

        let c0 = ring.add(
            &ring.add(&ring.mul(&self.pk.0, &u), &e1),
            &m.0.clone().scale(delta),
        );
        let c1 = ring.add(&ring.mul(&self.pk.1, &u), &e2);

        CipherText(c0, c1)
    }
}

pub struct Decrypt {
    params: BFVParams,
    sk: SecretKey,
}

impl Decrypt {
    pub fn decrypt(&self, c: &CipherText) -> Message {
        let ring = self.params.ring_q();
        let (q, t) = (self.params.Q as i64, self.params.t as i64);

        let mut v = ring.add(&c.0, &ring.mul(&c.1, &self.sk.0)); // reduced mod Q
        v.coef.iter_mut().for_each(|x| {
            if *x > q / 2 {
                *x -= q
            } // center lift
              // round(t·x / q) with integer arithmetic
            *x = (2 * t * *x + q).div_euclid(2 * q);
            *x = x.rem_euclid(t);
        });
        Message(v)
    }
}

#[derive(Debug, Clone)]
pub struct CipherText(Polynomial<i64>, Polynomial<i64>);

#[derive(Debug, Clone, PartialEq)]
pub struct Message(Polynomial<i64>);

pub fn generate_random_secret_key(len: usize) -> SecretKey {
    SecretKey::new(uniform_R2(len))
}

#[cfg(test)]
mod tests {
    use super::*;
    use tests::utils::*;

    /// The randomness is not seedable yet, so repeat with fresh keys instead.
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
    // not just `delta * m` plus a few units of noise.
    #[test]
    fn test_ciphertext_is_masked() {
        let q = test_params().Q as i64;
        let msg = Polynomial::from_slice(&[5, 4, 0, 0]);

        let (mut c0_masked, mut c1_masked) = (0, 0);
        for _ in 0..50 {
            let scheme = test_scheme();
            let delta = scheme.enc.params.delta();
            let c = scheme.enc.encrypt(Message(msg.clone()));

            let mask = c.0.clone() - msg.clone().scale(delta);
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
        let msg = Polynomial::from_slice(&[5, 4, 0, 0]);

        for _ in 0..TRIAL_COUNT {
            let (alice, mut bob) = (test_scheme(), test_scheme());
            // Ensure that we'll never hit the case where it's the same SK generated.
            while bob.dec.sk == alice.dec.sk {
                bob = test_scheme();
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

        for (t, q, a, b, expected) in cases {
            let scheme = scheme_with(t, q);
            let expected = Polynomial::from_slice(&expected);
            let ma = Message(Polynomial::from_slice(&a));
            let mb = Message(Polynomial::from_slice(&b));

            for _ in 0..TRIAL_COUNT {
                let ca = scheme.enc.encrypt(ma.clone());
                let cb = scheme.enc.encrypt(mb.clone());

                let product = scheme
                    .dec
                    .decrypt(&scheme.eval.mul(ca, cb, scheme.ek.clone()));

                assert_eq!(product.0, expected, "t = {t}, Q = {q}, {a:?} * {b:?}");
            }
        }
    }

    #[test]
    fn test_random_messages() {
        let (t, q) = (17usize, 525641usize);
        for _ in 0..TRIAL_COUNT {
            let scheme = scheme_with(t, q);
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
                "{a:?} + {b:?}"
            );

            let product = scheme.dec.decrypt(&scheme.eval.mul(ca, cb, scheme.ek));
            assert_eq!(product.0, plain_mul(&a, &b, 4, t as i64), "{a:?} * {b:?}");
        }
    }

    pub mod utils {
        use super::*;

        pub struct Scheme {
            pub enc: Encrypt,
            pub dec: Decrypt,
            pub eval: Evaluator,
            pub ek: EvalKey,
        }

        pub fn test_params() -> BFVParams {
            BFVParams {
                N: 4,
                t: 8,
                Q: 525641,
                P: 1024,
            }
        }

        pub fn scheme_with(t: usize, q: usize) -> Scheme {
            let params = BFVParams {
                t,
                Q: q,
                ..test_params()
            };

            let sk = generate_random_secret_key(params.N);
            let pk = PublicKey::generate_from_secret_key(&sk, params);
            let ek = EvalKey::generate_from_secret_key(&sk, params);

            Scheme {
                enc: Encrypt { params: params, pk },
                dec: Decrypt { params: params, sk },
                eval: Evaluator { params: params },
                ek,
            }
        }

        pub fn test_scheme() -> Scheme {
            let params = test_params();
            scheme_with(params.t, params.Q)
        }

        /// Reference product in Z_t[X]/(X^n + 1).
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

        /// Largest coefficient magnitude of `p` after centering mod `q`.
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
