use feanor_math;
use feanor_math::rings::multivariate::MultiVariatePolyRingEl;

struct QuotientPolynomialRing<T> {
    quotient: MultiVariatePolyRingEl<T>,
}

impl<T> QuotientPolynomialRing<T> {
    pub fn new(quotient: MultiVariatePolyRingEl<T>) -> Self {
        Self { quotient }
    }
}

// impl<T> RingBase for QuotientPolynomialRing<T> {}
//
// impl MultivariatePolyRing for QuotientPolynomialRing {}
