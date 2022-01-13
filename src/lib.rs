//! # black_scholes
//! A Black Scholes option pricing library.
use serde::Serialize;
use special::Error;
use std::f64::consts::{PI, SQRT_2};

fn cum_norm(x: f64) -> f64 {
    (x / SQRT_2).erf() * 0.5 + 0.5
}
fn inc_norm(x: f64) -> f64 {
    (-x.powi(2) / 2.0).exp() / (PI.sqrt() * SQRT_2)
}

fn d1(s: f64, k: f64, discount: f64, sqrt_maturity_sigma: f64) -> f64 {
    (s / (k * discount)).ln() / sqrt_maturity_sigma + 0.5 * sqrt_maturity_sigma
}
fn max_or_zero(v: f64) -> f64 {
    if v > 0.0 {
        v
    } else {
        0.0
    }
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
        let d1 = d1(s, k, discount, sqrt_maturity_sigma);
        s * cum_norm(d1) - k * discount * cum_norm(d1 - sqrt_maturity_sigma)
    } else {
        max_or_zero(s - k)
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
        let d1 = d1(s, k, discount, sqrt_maturity_sigma);
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
        let d1 = d1(s, k, discount, sqrt_maturity_sigma);
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
        let d1 = d1(s, k, discount, sqrt_maturity_sigma);
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
        let d1 = d1(s, k, discount, sqrt_maturity_sigma);
        -s * inc_norm(d1) * sigma / (2.0 * sqrt_t)
            - rate * k * discount * cum_norm(d1 - sqrt_maturity_sigma)
    } else {
        0.0
    }
}

/// Returns rho of a BS call option
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma=0.3;
/// let maturity=1.0;
/// let theta = black_scholes::call_rho(
///     stock, strike, rate, sigma, maturity
/// );
/// ```
pub fn call_rho(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    let sqrt_t = maturity.sqrt();
    let sqrt_maturity_sigma = sqrt_t * sigma;
    if sqrt_maturity_sigma > 0.0 {
        let discount = (-rate * maturity).exp();
        let d1 = d1(s, k, discount, sqrt_maturity_sigma);
        k * discount * maturity * cum_norm(d1 - sqrt_maturity_sigma)
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
        let d1 = d1(s, k, discount, sqrt_maturity_sigma);
        k * discount * cum_norm(sqrt_maturity_sigma - d1) - s * cum_norm(-d1)
    } else {
        max_or_zero(k - s)
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
        let d1 = d1(s, k, discount, sqrt_maturity_sigma);
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
        let d1 = d1(s, k, discount, sqrt_maturity_sigma);
        -s * inc_norm(d1) * sigma / (2.0 * sqrt_t)
            + rate * k * discount * cum_norm(-d1 + sqrt_maturity_sigma)
    } else {
        0.0
    }
}
/// Returns rho of a BS put option
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma=0.3;
/// let maturity=1.0;
/// let theta = black_scholes::put_rho(
///     stock, strike, rate, sigma, maturity
/// );
/// ```
pub fn put_rho(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    let sqrt_t = maturity.sqrt();
    let sqrt_maturity_sigma = sqrt_t * sigma;
    if sqrt_maturity_sigma > 0.0 {
        let discount = (-rate * maturity).exp();
        let d1 = d1(s, k, discount, sqrt_maturity_sigma);

        -1.0 * k * discount * maturity * cum_norm(-d1 + sqrt_maturity_sigma)
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

#[derive(Serialize)]
pub struct PricesAndGreeks {
    call_price: f64,
    call_delta: f64,
    call_gamma: f64,
    call_theta: f64,
    call_vega: f64,
    call_rho: f64,
    put_price: f64,
    put_delta: f64,
    put_gamma: f64,
    put_theta: f64,
    put_vega: f64,
    put_rho: f64,
}
/// Returns call and put prices and greeks.
/// Due to caching the complex computations
/// (such as N(d1)), this implementation is
/// faster if you need to obtain all the
/// information for a given stock price
/// and strike price.
///
/// # Examples
///
/// ```
/// let sigma = 0.3;
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let maturity = 1.0;
/// let all_prices_and_greeks = black_scholes::compute_all(
///     stock,
///     strike,
///     rate,
///     sigma,
///     maturity,
/// );
/// ```
pub fn compute_all(
    stock: f64,
    strike: f64,
    rate: f64,
    sigma: f64,
    maturity: f64,
) -> PricesAndGreeks {
    let discount = (-rate * maturity).exp();
    let sqrt_maturity = maturity.sqrt();
    let sqrt_maturity_sigma = sqrt_maturity * sigma;
    let k_discount = strike * discount;
    if sqrt_maturity_sigma > 0.0 {
        let d1 = d1(stock, strike, discount, sqrt_maturity_sigma);
        let d2 = d1 - sqrt_maturity_sigma;
        let cdf_d1 = cum_norm(d1);
        let cdf_d2 = cum_norm(d2);
        let pdf_d1 = inc_norm(d1);

        let call_price = stock * cdf_d1 - k_discount * cdf_d2;
        let call_delta = cdf_d1;
        let call_gamma = pdf_d1 / (stock * sqrt_maturity_sigma);
        let call_theta =
            -stock * pdf_d1 * sigma / (2.0 * sqrt_maturity) - rate * k_discount * cdf_d2;
        let call_vega = stock * pdf_d1 * sqrt_maturity_sigma / sigma;
        let call_rho = k_discount * maturity * cdf_d2;
        let put_price = call_price + k_discount - stock;
        let put_delta = cdf_d1 - 1.0;
        let put_gamma = call_gamma;
        let put_theta =
            -stock * pdf_d1 * sigma / (2.0 * sqrt_maturity) + rate * k_discount * (1.0 - cdf_d2);
        let put_vega = call_vega;
        let put_rho = -1.0 * k_discount * maturity * (1.0 - cdf_d2);
        PricesAndGreeks {
            call_price,
            call_delta,
            call_gamma,
            call_theta,
            call_vega,
            call_rho,
            put_price,
            put_delta,
            put_gamma,
            put_theta,
            put_vega,
            put_rho,
        }
    } else {
        PricesAndGreeks {
            call_price: max_or_zero(stock - strike),
            call_delta: if stock > strike { 1.0 } else { 0.0 },
            call_gamma: 0.0,
            call_theta: 0.0,
            call_vega: 0.0,
            call_rho: 0.0,
            put_price: max_or_zero(strike - stock),
            put_delta: if strike > stock { -1.0 } else { 0.0 },
            put_gamma: 0.0,
            put_theta: 0.0,
            put_vega: 0.0,
            put_rho: 0.0,
        }
    }
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

    macro_rules! assert_approx_eq {
        ($a:expr, $b:expr) => {{
            let (a, b) = (&$a, &$b);
            assert!(
                (*a - *b).abs() < 1.0e-6,
                "{} is not approximately equal to {}",
                *a,
                *b
            );
        }};
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
        assert_approx_eq!(call(5.0, 4.5, 0.05, 0.3, 1.0), 0.9848721043419868);
    }
    #[test]
    fn call_formula_works_close_maturity() {
        assert_approx_eq!(
            call(14341.0, 14000.0, 0.1, 0.2125, 0.25 / 365.0),
            341.95898982726794
        );
    }
    #[test]
    fn call_formula_works_with_zero_vol() {
        assert_eq!(call(5.0, 4.5, 0.05, 0.3, 0.0), 0.5);
    }
    #[test]
    fn put_formula_works() {
        assert_approx_eq!(put(5.0, 4.5, 0.05, 0.3, 1.0), 0.2654045145951993);
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
            let discount = (-rate * maturity).exp();
            let initial_guess = approximate_vol(price, s, k, rate, maturity);
            let cutoff = 0.000001;
            if price > cutoff && s - k * discount - price > cutoff {
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

    #[test]
    fn put_greeks_work() {
        let s = 550.88;
        let sigma = 0.37;
        let k = 510.0;
        let rate = 0.0;
        let maturity = 0.09;
        assert_approx_eq!(put_rho(s, k, rate, sigma, maturity), -11.996530249211213);
        assert_approx_eq!(put_theta(s, k, rate, sigma, maturity), -102.28760152696525);
        assert_approx_eq!(put_gamma(s, k, rate, sigma, maturity), 0.00492419827941365);
        assert_approx_eq!(put_vega(s, k, rate, sigma, maturity), 49.761535877983086);
        assert_approx_eq!(put_delta(s, k, rate, sigma, maturity), -0.22658184828282102);
    }

    #[test]
    fn call_greeks_works() {
        let s = 550.88;
        let sigma = 0.37;
        let k = 510.0;
        let rate = 0.0;
        let maturity = 0.09;
        assert_approx_eq!(call_rho(s, k, rate, sigma, maturity), 33.90346975078879);
        assert_approx_eq!(call_theta(s, k, rate, sigma, maturity), -102.28760152696525);
        assert_approx_eq!(call_gamma(s, k, rate, sigma, maturity), 0.00492419827941365);
        assert_approx_eq!(call_vega(s, k, rate, sigma, maturity), 49.761535877983086);
        assert_approx_eq!(call_delta(s, k, rate, sigma, maturity), 0.773418151717179);
    }
    #[test]
    fn compute_all_works() {
        let s = 550.88;
        let sigma = 0.37;
        let k = 510.0;
        let rate = 0.0;
        let maturity = 0.09;
        let PricesAndGreeks {
            call_price,
            call_delta,
            call_gamma,
            call_theta,
            call_vega,
            call_rho,
            put_price,
            put_delta,
            put_gamma,
            put_theta,
            put_vega,
            put_rho,
        } = compute_all(s, k, rate, sigma, maturity);
        assert_approx_eq!(call_price, call(s, k, rate, sigma, maturity));
        assert_approx_eq!(call_delta, 0.773418151717179);
        assert_approx_eq!(call_gamma, 0.00492419827941365);
        assert_approx_eq!(call_theta, -102.28760152696525);
        assert_approx_eq!(call_vega, 49.761535877983086);
        assert_approx_eq!(call_rho, 33.90346975078879);

        assert_approx_eq!(put_price, put(s, k, rate, sigma, maturity));
        assert_approx_eq!(put_delta, -0.22658184828282102);
        assert_approx_eq!(put_gamma, 0.00492419827941365);
        assert_approx_eq!(put_theta, -102.28760152696525);
        assert_approx_eq!(put_vega, 49.761535877983086);
        assert_approx_eq!(put_rho, -11.996530249211213);
    }

    #[test]
    fn compute_all_works_rate() {
        let s = 550.88;
        let sigma = 0.37;
        let k = 510.0;
        let rate = 0.05;
        let maturity = 0.09;
        let result = compute_all(s, k, rate, sigma, maturity);
        assert_approx_eq!(result.call_price, call(s, k, rate, sigma, maturity));
        assert_approx_eq!(result.call_delta, call_delta(s, k, rate, sigma, maturity));
        assert_approx_eq!(result.call_gamma, call_gamma(s, k, rate, sigma, maturity));
        assert_approx_eq!(result.call_theta, call_theta(s, k, rate, sigma, maturity));
        assert_approx_eq!(result.call_vega, call_vega(s, k, rate, sigma, maturity));
        assert_approx_eq!(result.call_rho, call_rho(s, k, rate, sigma, maturity));

        assert_approx_eq!(result.put_price, put(s, k, rate, sigma, maturity));
        assert_approx_eq!(result.put_delta, put_delta(s, k, rate, sigma, maturity));
        assert_approx_eq!(result.put_gamma, put_gamma(s, k, rate, sigma, maturity));
        assert_approx_eq!(result.put_theta, put_theta(s, k, rate, sigma, maturity));
        assert_approx_eq!(result.put_vega, put_vega(s, k, rate, sigma, maturity));
        assert_approx_eq!(result.put_rho, put_rho(s, k, rate, sigma, maturity));
    }
}
