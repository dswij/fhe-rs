use criterion::{black_box, criterion_group, criterion_main, Criterion};

use fhe_rs::dft::fft;
use num::complex::Complex;

fn criterion_benchmark(c: &mut Criterion) {
    let input = [Complex::new(1.0, 0.0); 100000];
    c.bench_function("fft", |b| b.iter(|| fft(&input)));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
