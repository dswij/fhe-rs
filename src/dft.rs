use num::complex::Complex;

use std::f64::consts::PI;

/// Implementation of Cooley-Tukey FFT
pub fn fft(w: &[Complex<f64>]) -> Vec<Complex<f64>> {
    if w.len() == 1 {
        return w.iter().cloned().collect();
    }

    let odd: Vec<_> = w.iter().step_by(2).cloned().collect();
    let even: Vec<_> = w[1..].iter().step_by(2).cloned().collect();

    let F_odd = fft(&odd);
    let F_even = fft(&even);

    let mut result: Vec<Complex<f64>> = vec![Complex::default(); w.len()];

    for j in 0..w.len() / 2 {
        let w_n = Complex::from_polar(1.0, (j as f64) * -2.0 * PI / (w.len() as f64));
        result[j] = F_even[j] + w_n * F_odd[j];
        result[j + w.len() / 2] = F_even[j] - w_n * F_odd[j];
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fft_all_ones() {
        let input = [Complex::new(1.0, 0.0); 8];

        let mut expected = [Complex::new(0.0, 0.0); 8];
        expected[0] = Complex::new(8.0, 0.0);

        assert_eq!(fft(&input), expected);
    }

    #[test]
    fn test_fft_all_zeroes() {}

    #[test]
    fn test_fft_uniform() {}
}
