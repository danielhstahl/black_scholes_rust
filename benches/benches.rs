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
