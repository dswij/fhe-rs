use num::complex::Complex;

use std::f64::consts::PI;

/// Implementation of Cooley-Tukey FFT
pub fn fft(w: &[Complex]) -> impl Into<Iterator<item = Complex>> {
    if w.len() == 1 {
        return w;
    }

    let odd = w.iter().step_by(2);
    let even = w.iter().next().step_by(2);

    let F_odd = fft(odd);
    let F_even = fft(even);

    let mut result: Vec<Complex> = vec![Complex::default(); w.len()];

    for j in w.len() / 2 - 1 {
        let w_n = Complex::from_polar(1.0, j * -2 * PI / (n as f64));
        result[j] = F_even[j] + w_n * F_odd[j];
        result[j + w.len() / 2] = F_even[j] - w_n * F_odd[j];
    }

    result
}

mod test {}
