#![feature(test)]
extern crate test;
use test::Bencher;
#[bench]
fn bench_call_price(b: &mut Bencher) {
    let r = 0.05;
    let sig = 0.3;
    let t = 1.0;
    let asset = 50.0;
    let k = 50.0;
    b.iter(|| black_scholes::call(asset, k, r, sig, t))
}

#[bench]
fn bench_all_price(b: &mut Bencher) {
    let r = 0.05;
    let sig = 0.3;
    let t = 1.0;
    let asset = 50.0;
    let k = 50.0;
    b.iter(|| black_scholes::compute_all(asset, k, r, sig, t))
}

#[bench]
fn bench_all_price_no_cache(b: &mut Bencher) {
    let r = 0.05;
    let sig = 0.3;
    let t = 1.0;
    let asset = 50.0;
    let k = 50.0;
    b.iter(|| {
        black_scholes::call(asset, k, r, sig, t);
        black_scholes::call_delta(asset, k, r, sig, t);
        black_scholes::call_gamma(asset, k, r, sig, t);
        black_scholes::call_vega(asset, k, r, sig, t);
        black_scholes::call_theta(asset, k, r, sig, t);
        black_scholes::call_rho(asset, k, r, sig, t);

        black_scholes::put(asset, k, r, sig, t);
        black_scholes::put_delta(asset, k, r, sig, t);
        black_scholes::put_gamma(asset, k, r, sig, t);
        black_scholes::put_vega(asset, k, r, sig, t);
        black_scholes::put_theta(asset, k, r, sig, t);
        black_scholes::put_rho(asset, k, r, sig, t);
    })
}
