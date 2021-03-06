//! # black_scholes
//! A Black Scholes option pricing library.

use special::Error;
use std::f64::consts::{PI, SQRT_2};

fn cum_norm(x: f64) -> f64 {
    (x / SQRT_2).erf() * 0.5 + 0.5
}
fn inc_norm(x: f64) -> f64 {
    (-x.powi(2) / 2.0).exp() / (PI.sqrt() * SQRT_2)
}

/// Returns BS call option formula with discount and volatility already computed.
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let discount = 0.99;
/// let sigma = 0.3;
/// let maturity:f64 = 2.0;
/// let sqrt_maturity_sigma = sigma*maturity.sqrt();
/// let price = black_scholes::call_discount(
///     stock, strike, discount,
///     sqrt_maturity_sigma
/// );
/// ```
pub fn call_discount(s: f64, k: f64, discount: f64, sqrt_maturity_sigma: f64) -> f64 {
    if sqrt_maturity_sigma > 0.0 {
        let d1 = (s / (k * discount)).ln() / sqrt_maturity_sigma + 0.5 * sqrt_maturity_sigma;
        s * cum_norm(d1) - k * discount * cum_norm(d1 - sqrt_maturity_sigma)
    } else {
        if s > k {
            s - k
        } else {
            0.0
        }
    }
}

/// Returns standard BS call option formula.
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma=0.3;
/// let maturity=1.0;
/// assert_eq!(0.9848721043419868, black_scholes::call(stock, strike, rate, sigma, maturity));
/// ```
pub fn call(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    call_discount(s, k, (-rate * maturity).exp(), maturity.sqrt() * sigma)
}

/// Returns delta of a BS call option
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma=0.3;
/// let maturity=1.0;
/// let delta = black_scholes::call_delta(
///     stock, strike, rate, sigma, maturity
/// );
/// ```
pub fn call_delta(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    let sqrt_maturity_sigma = maturity.sqrt() * sigma;
    if sqrt_maturity_sigma > 0.0 {
        let discount = (-rate * maturity).exp();
        let d1 = (s / (k * discount)).ln() / sqrt_maturity_sigma + 0.5 * sqrt_maturity_sigma;
        cum_norm(d1)
    } else {
        if s > k {
            1.0
        } else {
            0.0
        }
    }
}

/// Returns gamma of a BS call option
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma=0.3;
/// let maturity=1.0;
/// let gamma = black_scholes::call_gamma(
///     stock, strike, rate, sigma, maturity
/// );
/// ```
pub fn call_gamma(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    let sqrt_maturity_sigma = maturity.sqrt() * sigma;
    if sqrt_maturity_sigma > 0.0 {
        let discount = (-rate * maturity).exp();
        let d1 = (s / (k * discount)).ln() / sqrt_maturity_sigma + 0.5 * sqrt_maturity_sigma;
        inc_norm(d1) / (s * sqrt_maturity_sigma)
    } else {
        0.0
    }
}
/// Returns vega of a BS call option
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma=0.3;
/// let maturity=1.0;
/// let vega = black_scholes::call_vega(
///     stock, strike, rate, sigma, maturity
/// );
/// ```
pub fn call_vega(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    let sqrt_maturity_sigma = maturity.sqrt() * sigma;
    if sqrt_maturity_sigma > 0.0 {
        let discount = (-rate * maturity).exp();
        let d1 = (s / (k * discount)).ln() / sqrt_maturity_sigma + 0.5 * sqrt_maturity_sigma;
        s * inc_norm(d1) * sqrt_maturity_sigma / sigma
    } else {
        0.0
    }
}
/// Returns theta of a BS call option
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma=0.3;
/// let maturity=1.0;
/// let theta = black_scholes::call_theta(
///     stock, strike, rate, sigma, maturity
/// );
/// ```
pub fn call_theta(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    let sqrt_t = maturity.sqrt();
    let sqrt_maturity_sigma = sqrt_t * sigma;
    if sqrt_maturity_sigma > 0.0 {
        let discount = (-rate * maturity).exp();
        let d1 = (s / (k * discount)).ln() / sqrt_maturity_sigma + 0.5 * sqrt_maturity_sigma;
        -s * inc_norm(d1) * sigma / (2.0 * sqrt_t)
            - rate * k * discount * cum_norm(d1 - sqrt_maturity_sigma)
    } else {
        0.0
    }
}

/// Returns BS put option formula with discount and volatility already computed.
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let discount = 0.99;
/// let sigma = 0.3;
/// let maturity:f64 = 2.0;
/// let sqrt_maturity_sigma = sigma*maturity.sqrt();
/// let price = black_scholes::put_discount(
///     stock, strike, discount,
///     sqrt_maturity_sigma
/// );
/// ```
pub fn put_discount(s: f64, k: f64, discount: f64, sqrt_maturity_sigma: f64) -> f64 {
    if sqrt_maturity_sigma > 0.0 {
        let d1 = (s / (k * discount)).ln() / sqrt_maturity_sigma + 0.5 * sqrt_maturity_sigma;
        k * discount * cum_norm(sqrt_maturity_sigma - d1) - s * cum_norm(-d1)
    } else {
        if k > s {
            k - s
        } else {
            0.0
        }
    }
}

/// Returns BS put option formula.
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma = 0.3;
/// let maturity = 1.0;
/// assert_eq!(0.2654045145951993, black_scholes::put(stock, strike, rate, sigma, maturity));
/// ```
pub fn put(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    put_discount(s, k, (-rate * maturity).exp(), maturity.sqrt() * sigma)
}

/// Returns delta of a BS put option
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma=0.3;
/// let maturity=1.0;
/// let delta = black_scholes::put_delta(
///     stock, strike, rate, sigma, maturity
/// );
/// ```
pub fn put_delta(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    let sqrt_maturity_sigma = maturity.sqrt() * sigma;
    if sqrt_maturity_sigma > 0.0 {
        let discount = (-rate * maturity).exp();
        let d1 = (s / (k * discount)).ln() / sqrt_maturity_sigma + 0.5 * sqrt_maturity_sigma;
        return cum_norm(d1) - 1.0;
    } else {
        return if k > s { -1.0 } else { 0.0 };
    }
}
/// Returns gamma of a BS put option
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma=0.3;
/// let maturity=1.0;
/// let gamma = black_scholes::put_gamma(
///     stock, strike, rate, sigma, maturity
/// );
/// ```
pub fn put_gamma(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    call_gamma(s, k, rate, sigma, maturity) //same as call
}

/// Returns vega of a BS put option
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma=0.3;
/// let maturity=1.0;
/// let vega = black_scholes::put_vega(
///     stock, strike, rate, sigma, maturity
/// );
/// ```
pub fn put_vega(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    call_vega(s, k, rate, sigma, maturity) //same as call
}

/// Returns theta of a BS put option
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma=0.3;
/// let maturity=1.0;
/// let theta = black_scholes::put_theta(
///     stock, strike, rate, sigma, maturity
/// );
/// ```
pub fn put_theta(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    let sqrt_t = maturity.sqrt();
    let sqrt_maturity_sigma = sqrt_t * sigma;
    if sqrt_maturity_sigma > 0.0 {
        let discount = (-rate * maturity).exp();
        let d1 = (s / (k * discount)).ln() / sqrt_maturity_sigma + 0.5 * sqrt_maturity_sigma;
        -s * inc_norm(d1) * sigma / (2.0 * sqrt_t)
            + rate * k * discount * cum_norm(-d1 + sqrt_maturity_sigma)
    } else {
        0.0
    }
}
const SQRT_TWO_PI: f64 = 2.0 * std::f64::consts::SQRT_2 / std::f64::consts::FRAC_2_SQRT_PI;
//Corrado and Miller (1996)
fn approximate_vol(price: f64, s: f64, k: f64, rate: f64, maturity: f64) -> f64 {
    let discount = (-rate * maturity).exp();
    let x = k * discount;
    let coef = SQRT_TWO_PI / (s + x);
    let helper_1 = s - x;
    let c1 = price - helper_1 * 0.5;
    let c2 = c1.powi(2);
    let c3 = helper_1.powi(2) / std::f64::consts::PI;
    let bridge_1 = c2 - c3;
    let bridge_m = if bridge_1 > 0.0 { bridge_1.sqrt() } else { 0.0 };
    coef * (c1 + bridge_m) / maturity.sqrt()
}
/// Returns implied volatility from a call option with initial guess
///
/// # Examples
///
/// ```
/// let price = 1.0;
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let maturity = 1.0;
/// let initial_guess = 0.3;
/// let iv = black_scholes::call_iv_guess(
///     price, stock, strike, rate,
///     maturity, initial_guess
/// ).unwrap();
/// ```
pub fn call_iv_guess(
    price: f64,
    s: f64,
    k: f64,
    rate: f64,
    maturity: f64,
    initial_guess: f64,
) -> Result<f64, f64> {
    let obj_fn = |sigma| call(s, k, rate, sigma, maturity) - price;
    let dfn = |sigma| call_vega(s, k, rate, sigma, maturity);
    let precision = 0.000001;
    let iterations = 10000;
    nrfind::find_root(&obj_fn, &dfn, initial_guess, precision, iterations)
}
/// Returns implied volatility from a call option
///
/// # Examples
///
/// ```
/// let price = 1.0;
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let maturity = 1.0;
/// let iv = black_scholes::call_iv(
///     price, stock, strike, rate,
///     maturity
/// ).unwrap();
/// ```
pub fn call_iv(price: f64, s: f64, k: f64, rate: f64, maturity: f64) -> Result<f64, f64> {
    let initial_guess = approximate_vol(price, s, k, rate, maturity);
    call_iv_guess(price, s, k, rate, maturity, initial_guess)
}

/// Returns implied volatility from a put option with initial guess
///
/// # Examples
///
/// ```
/// let price = 0.3;
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let maturity = 1.0;
/// let initial_guess = 0.3;
/// let iv = black_scholes::put_iv_guess(
///     price, stock, strike, rate,
///     maturity, initial_guess
/// ).unwrap();
/// ```
pub fn put_iv_guess(
    price: f64,
    s: f64,
    k: f64,
    rate: f64,
    maturity: f64,
    initial_guess: f64,
) -> Result<f64, f64> {
    let obj_fn = |sigma| put(s, k, rate, sigma, maturity) - price;
    let dfn = |sigma| put_vega(s, k, rate, sigma, maturity);
    let precision = 0.000001;
    let iterations = 10000;
    nrfind::find_root(&obj_fn, &dfn, initial_guess, precision, iterations)
}
/// Returns implied volatility from a put option
///
/// # Examples
///
/// ```
/// let price = 0.3;
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let maturity = 1.0;
/// let initial_guess = 0.3;
/// let iv = black_scholes::put_iv(
///     price, stock, strike, rate,
///     maturity
/// ).unwrap();
/// ```
pub fn put_iv(price: f64, s: f64, k: f64, rate: f64, maturity: f64) -> Result<f64, f64> {
    let c_price = price + s - k * (-rate * maturity).exp();
    let initial_guess = approximate_vol(c_price, s, k, rate, maturity);
    put_iv_guess(price, s, k, rate, maturity, initial_guess)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::*;
    use rand::distributions::{Distribution, Uniform};
    use rand::{SeedableRng, StdRng};
    fn get_rng_seed(seed: [u8; 32]) -> StdRng {
        SeedableRng::from_seed(seed)
    }
    fn get_over_region(lower: f64, upper: f64, rand: f64) -> f64 {
        lower + (upper - lower) * rand
    }
    #[test]
    fn sqrt_two_pi_is_right() {
        assert_abs_diff_eq!(
            SQRT_TWO_PI,
            (2.0 * std::f64::consts::PI).sqrt(),
            epsilon = 0.000000001
        );
    }
    #[test]
    fn call_formula_works() {
        assert_eq!(call(5.0, 4.5, 0.05, 0.3, 1.0), 0.9848721043419868);
    }
    #[test]
    fn call_formula_works_with_zero_vol() {
        assert_eq!(call(5.0, 4.5, 0.05, 0.3, 0.0), 0.5);
    }
    #[test]
    fn put_formula_works() {
        assert_eq!(put(5.0, 4.5, 0.05, 0.3, 1.0), 0.2654045145951993);
    }
    #[test]
    fn put_formula_works_with_zero_vol() {
        assert_eq!(put(5.0, 4.5, 0.05, 0.3, 0.0), 0.0);
    }
    #[test]
    fn call_iv_works() {
        let sigma = 0.2;
        let initial_guess = 0.5;
        let s = 5.0;
        let k = 4.5;
        let rate = 0.05;
        let maturity = 1.0;
        let price = call(s, k, rate, sigma, maturity);
        assert_abs_diff_eq!(
            call_iv_guess(price, s, k, rate, maturity, initial_guess).unwrap(),
            sigma,
            epsilon = 0.00000001
        );
    }
    #[test]
    fn call_iv_approx() {
        let sigma = 0.2;
        let s = 5.0;
        let k = 4.5;
        let rate = 0.05;
        let maturity = 1.0;
        let price = call(s, k, rate, sigma, maturity);
        let approx_vol = approximate_vol(price, s, k, rate, maturity);
        println!("approx: {}", approx_vol);
        assert_abs_diff_eq!(sigma, approx_vol, epsilon = 0.01);
    }
    #[test]
    fn call_iv_works_with_broad_set_of_numbers() {
        let seed: [u8; 32] = [2; 32];
        let mut rng_seed = get_rng_seed(seed);
        let uniform = Uniform::new(0.0f64, 1.0);
        let num_total: usize = 10000;

        (0..num_total).for_each(|_| {
            let s = 1.0;
            let k = get_over_region(0.3, 3.0, uniform.sample(&mut rng_seed));
            let sigma = get_over_region(0.1, 2.0, uniform.sample(&mut rng_seed));
            let rate = 0.0247;
            let maturity = 0.7599;
            let price = call(s, k, rate, sigma, maturity);
            let initial_guess = approximate_vol(price, s, k, rate, maturity);
            //println!("s: {}, k: {}, sigma: {}, price: {}, initial_guess: {}", s, k, sigma, price, initial_guess);
            if price > 0.000001 {
                let _iv = call_iv_guess(price, s, k, rate, maturity, initial_guess).unwrap();
            }
        })
    }
    #[test]
    fn call_iv_works_with_difficult() {
        let s = 0.43065239380643594;
        let k = 0.5016203266170813;
        let sigma = 0.4192621453186373;
        let rate = 0.0247;
        let maturity = 0.7599;
        let price = call(s, k, rate, sigma, maturity);
        println!("s: {}, k: {}, sigma: {}, price: {}", s, k, sigma, price);
        let _iv = call_iv(price, s, k, rate, maturity).unwrap();
    }
    #[test]
    fn put_iv_works() {
        let sigma = 0.2;
        //let initial_guess=0.195;
        let s = 5.0;
        let k = 4.5;
        let rate = 0.05;
        let maturity = 1.0;
        let price = put(s, k, rate, sigma, maturity);
        assert_abs_diff_eq!(
            put_iv(price, s, k, rate, maturity).unwrap(),
            sigma,
            epsilon = 0.00000001
        );
    }
    #[test]
    fn call_iv_returns_err_if_no_possible_solution() {
        let price = 50.275;
        let s = 274.525;
        let k = 225.000;
        let rate = 0.0244;
        let maturity = 0.156;
        assert!(call_iv(price, s, k, rate, maturity).is_err());
    }
}
