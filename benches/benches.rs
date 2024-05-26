#![feature(test)]
extern crate test;
use test::{black_box, Bencher};

/// Convenience alias for [`black_box`].
const BB: fn(f64) -> f64 = black_box;

#[bench]
fn bench_call_price(b: &mut Bencher) {
    let r = 0.05;
    let sig = 0.3;
    let t = 1.0;
    let asset = 50.0;
    let k = 50.0;
    b.iter(|| black_scholes::call(BB(asset), BB(k), BB(r), BB(sig), BB(t)))
}

#[bench]
fn bench_all_price(b: &mut Bencher) {
    let r = 0.05;
    let sig = 0.3;
    let t = 1.0;
    let asset = 50.0;
    let k = 50.0;
    b.iter(|| black_scholes::compute_all(BB(asset), BB(k), BB(r), BB(sig), BB(t)))
}

#[bench]
fn bench_all_price_no_cache(b: &mut Bencher) {
    let r = 0.05;
    let sig = 0.3;
    let t = 1.0;
    let asset = 50.0;
    let k = 50.0;
    b.iter(|| {
        black_scholes::call(BB(asset), BB(k), BB(r), BB(sig), BB(t));
        black_scholes::call_delta(BB(asset), BB(k), BB(r), BB(sig), BB(t));
        black_scholes::call_gamma(BB(asset), BB(k), BB(r), BB(sig), BB(t));
        black_scholes::call_vega(BB(asset), BB(k), BB(r), BB(sig), BB(t));
        black_scholes::call_theta(BB(asset), BB(k), BB(r), BB(sig), BB(t));
        black_scholes::call_rho(BB(asset), BB(k), BB(r), BB(sig), BB(t));

        black_scholes::put(BB(asset), BB(k), BB(r), BB(sig), BB(t));
        black_scholes::put_delta(BB(asset), BB(k), BB(r), BB(sig), BB(t));
        black_scholes::put_gamma(BB(asset), BB(k), BB(r), BB(sig), BB(t));
        black_scholes::put_vega(BB(asset), BB(k), BB(r), BB(sig), BB(t));
        black_scholes::put_theta(BB(asset), BB(k), BB(r), BB(sig), BB(t));
        black_scholes::put_rho(BB(asset), BB(k), BB(r), BB(sig), BB(t));
    })
}
