use crate::Polynomial;

#[derive(Debug)]
pub struct QuotientRing {
    pub quotient: Polynomial<i64>,
}

impl QuotientRing {
    pub fn add(&self, lhs: &Polynomial<i64>, rhs: &Polynomial<i64>) -> Polynomial<i64> {
        lhs + rhs
    }

    pub fn mul(&self, lhs: &Polynomial<i64>, rhs: &Polynomial<i64>) -> Polynomial<i64> {
        let (_, rem) = (lhs * rhs).div_rem(self.quotient.clone());

        rem
    }
}

#[derive(Debug)]
pub struct ResidualQuotientRing {
    pub ring: QuotientRing,
    pub bound: usize,
    t: i64,
}

impl ResidualQuotientRing {
    pub fn cyclotomic(degree: usize, bound: usize) -> Self {
        let mut quotient = vec![0; degree + 1];
        quotient[0] = 1;
        quotient[degree] = 1;
        let ring = QuotientRing {
            quotient: Polynomial::new(quotient),
        };

        Self {
            ring,
            bound,
            t: bound.div_ceil(2).try_into().unwrap(),
        }
    }

    pub fn add(&self, lhs: &Polynomial<i64>, rhs: &Polynomial<i64>) -> Polynomial<i64> {
        let mut result = self.ring.add(lhs, rhs);
        let bound = self.bound.try_into().unwrap();
        result.coef_mod(bound);

        result
    }

    pub fn mul(&self, lhs: &Polynomial<i64>, rhs: &Polynomial<i64>) -> Polynomial<i64> {
        let mut result = self.ring.mul(lhs, rhs);
        let bound = self.bound.try_into().unwrap();
        result.coef_mod(bound);

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quotient_ring_mul() {
        let p1 = Polynomial::from_slice(&[160, 45, 34, 0]);
        let p2 = Polynomial::from_slice(&[160, 45, 34, 0]);
        let ring = QuotientRing {
            quotient: Polynomial::from_slice(&[1, 0, 0, 0, 1]),
        };

        let result = ring.mul(&p1, &p2);

        let expected = Polynomial::from_slice(&[24444, 14400, 12905, 3060]);
        assert_eq!(result, expected);
    }
}
