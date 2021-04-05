#[macro_use]
extern crate criterion;
extern crate black_scholes;
use criterion::Criterion;

fn bench_call_price(c: &mut Criterion) {
    let r = 0.05;
    let sig = 0.3;
    let t = 1.0;
    let asset = 50.0;
    let k = 50.0;
    c.bench_function("call price", move |b| {
        b.iter(|| black_scholes::call(asset, k, r, sig, t))
    });
}

criterion_group!(benches, bench_call_price);
criterion_main!(benches);
