//! Modulus chain module that is required for modulus switching.
//!
//! This module covers the generation of such chain.

// {q_0, …, q_L}
#[derive(Debug, Clone, PartialEq)]
pub struct ModulusChain {
    primes: Vec<usize>,
}

impl ModulusChain {
    pub fn new(primes: Vec<usize>) -> Self {
        assert!(!primes.is_empty(), "primes can't be empty");

        // Check coprimality
        for (i, a) in primes.iter().enumerate() {
            for b in &primes[i + 1..] {
                assert_eq!(gcd(*a, *b), 1, "moduli {a} and {b} are not coprime");
            }
        }

        // Check that we won't overflow
        primes
            .iter()
            .try_fold(1usize, |acc, p| acc.checked_mul(*p))
            .expect("the product of the moduli overflows usize");

        Self::new_unchecked(primes)
    }

    pub fn new_unchecked(primes: Vec<usize>) -> Self {
        Self { primes }
    }

    // Primes `q_i ≡ 1 (mod lcm(t, 2n))` of about `bits[i]` bits each.
    pub fn generate(t: usize, n: usize, bits: &[u32]) -> Self {
        let m = lcm(t, 2 * n);
        let mut primes = Vec::with_capacity(bits.len());
        for b in bits {
            primes.push(find_prime(*b, m, &primes));
        }

        Self::new(primes)
    }

    // The highest level, i.e. the number of `mod_switch`es available.
    pub fn top_level(&self) -> usize {
        self.primes.len() - 1
    }

    // `Q_level = q_0 ⋯ q_level`.
    pub fn modulus(&self, level: usize) -> usize {
        self.primes[..=level].iter().product()
    }

    // The fresh-ciphertext modulus `Q_L`.
    pub fn top(&self) -> usize {
        self.modulus(self.top_level())
    }

    // Get prime at level.
    pub fn prime(&self, level: usize) -> usize {
        self.primes[level]
    }
}

// The largest prime `k·m + 1 ≤ 2^bits + 1` that is not in `skip`.
// TODO: find a better way to generate primes. Must be something out there
pub fn find_prime(bits: u32, m: usize, skip: &[usize]) -> usize {
    let mut k = (1usize << bits) / m;
    while k > 0 {
        let candidate = k * m + 1;
        if is_prime(candidate as u64) && !skip.contains(&candidate) {
            return candidate;
        }
        k -= 1;
    }

    panic!("no prime = 1 (mod {m}) of {bits} bits");
}

// Miller–Rabin
pub fn is_prime(n: u64) -> bool {
    const BASES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

    if n < 2 {
        return false;
    }
    for p in BASES {
        if n % p == 0 {
            return n == p;
        }
    }

    let (mut d, mut s) = (n - 1, 0);
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }

    'bases: for a in BASES {
        let mut x = pow_mod(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 0..s - 1 {
            x = mul_mod(x, x, n);
            if x == n - 1 {
                continue 'bases;
            }
        }
        return false;
    }

    true
}

// Helper to prevent overflowing
// TODO: read on how to handle 128-bits primitives in Rust.
fn mul_mod(a: u64, b: u64, m: u64) -> u64 {
    (a as u128 * b as u128 % m as u128) as u64
}

// Same as mul_mod, but for POW!
fn pow_mod(mut base: u64, mut exp: u64, m: u64) -> u64 {
    let mut result = 1;
    base %= m;
    while exp > 0 {
        if exp & 1 == 1 {
            result = mul_mod(result, base, m);
        }
        base = mul_mod(base, base, m);
        exp >>= 1;
    }

    result
}

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

fn lcm(a: usize, b: usize) -> usize {
    a / gcd(a, b) * b
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_prime() {
        let primes = [2, 3, 5, 457, 1009, 525641, 4_294_967_291];
        let composites = [
            0,
            1,
            4,
            1025,
            1089,
            561, /* Carmichael */
            4_294_967_297,
        ];

        assert!(primes.iter().all(|p| is_prime(*p)));
        assert!(composites.iter().all(|c| !is_prime(*c)));
    }

    #[test]
    fn test_generate_matches_hand_search() {
        // m = lcm(8, 8) = 8. Walking down from 2^10 and 2^9 gives 1009 and 457
        // (1025, 1017, 513, 505, … are composite).
        let chain = ModulusChain::generate(8, 4, &[10, 9]);
        assert_eq!(chain.primes, vec![1009, 457]);

        // m = lcm(17, 8) = 136.
        let chain = ModulusChain::generate(17, 4, &[10, 9]);
        assert_eq!(chain.primes, vec![953, 409]);
    }

    #[test]
    fn test_generate_constraints() {
        for (t, n) in [(8, 4), (17, 4), (8, 16)] {
            let chain = ModulusChain::generate(t, n, &[12, 10, 10, 10]);

            assert_eq!(chain.top_level(), 3);
            for level in 0..=chain.top_level() {
                let q = chain.prime(level);
                assert!(is_prime(q as u64));
                assert_eq!(q % t, 1, "q = {q} must be 1 mod t = {t}");
                assert_eq!(q % (2 * n), 1, "q = {q} must be 1 mod 2n = {}", 2 * n);
            }
        }
    }

    #[test]
    fn test_generate_primes_are_distinct() {
        // Three requests for the same size must not return the same prime.
        let chain = ModulusChain::generate(8, 4, &[9, 9, 9]);
        let mut primes = chain.primes.clone();
        primes.sort();
        primes.dedup();

        assert_eq!(primes.len(), 3);
    }

    #[test]
    fn test_modulus_is_prefix_product() {
        let chain = ModulusChain::new(vec![3, 5, 7]);

        assert_eq!(chain.modulus(0), 3);
        assert_eq!(chain.modulus(1), 15);
        assert_eq!(chain.modulus(2), 105);
        assert_eq!(chain.top(), 105);
        assert_eq!(chain.top() % chain.modulus(1), 0);
    }

    #[test]
    #[should_panic(expected = "not coprime")]
    fn test_new_rejects_non_coprime() {
        ModulusChain::new(vec![17, 34]);
    }
}
