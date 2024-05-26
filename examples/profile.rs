use std::hint::black_box;

/// Convenience alias for [`black_box`].
const BB: fn(f64) -> f64 = black_box;

fn main() {
    let r = 0.05;
    let sig = 0.3;
    let t = 1.0;
    let asset = 50.0;
    let k = 50.0;

    for _ in 0..100_000 {
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
    }
}
