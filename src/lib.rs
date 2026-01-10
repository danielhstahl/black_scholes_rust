//! # black_scholes
//! A Black Scholes option pricing library.
//!
//! This library provides implementations of the Black-Scholes option pricing model
//! along with the Black 76 model for options on futures and the Black-Scholes-Merton
//! model for options on dividend-paying stocks.
//!
//! ## Features
//! - European call and put option pricing
//! - Full set of Greeks (Delta, Gamma, Theta, Vega, Rho, Vanna, Vomma, Charm)
//! - Implied volatility calculation
//! - Support for multiple option models
//! - High-performance calculations with caching
use serde::Serialize;
use special::Error;
use std::f64::consts::{FRAC_1_PI, FRAC_1_SQRT_2, FRAC_2_SQRT_PI, SQRT_2};

/// 1/sqrt(2π)
#[allow(clippy::excessive_precision)]
const FRAC_1_SQRT_2PI: f64 = 0.3989422804014326779399460599343818684758586311649346576659258296;

/// Cumulative distribution function of standard normal distribution.
///
/// # Arguments
/// * `x` - Input value
///
/// # Returns
/// The cumulative probability up to `x`
fn cum_norm(x: f64) -> f64 {
    (x * FRAC_1_SQRT_2).error() * 0.5 + 0.5
}

/// Probability density function of standard normal distribution.
///
/// # Arguments
/// * `x` - Input value
///
/// # Returns
/// The probability density at `x`
fn inc_norm(x: f64) -> f64 {
    (-x.powi(2) * 0.5).exp() * FRAC_1_SQRT_2PI
}

/// Helper function to calculate d1 in the Black-Scholes formula.
///
/// # Arguments
/// * `s` - Stock price
/// * `k` - Strike price
/// * `discount` - Discount factor (e^(-rate * maturity))
/// * `sqrt_maturity_sigma` - Square root of (maturity * sigma^2)
///
/// # Returns
/// The calculated d1 value
fn d1(s: f64, k: f64, discount: f64, sqrt_maturity_sigma: f64) -> f64 {
    // equiv. to : ((s / k).ln() + (rate + 0.5 * sigma.powi(2)) * maturity) / sqrt_maturity_sigma
    (s / (k * discount)).ln() / sqrt_maturity_sigma + 0.5 * sqrt_maturity_sigma
}

#[inline(always)]
fn max_or_zero(v: f64) -> f64 {
    v.max(0.0)
}

/// Returns BS call option formula with discount and volatility already computed.
///
/// # Arguments
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `discount` - Discount factor (e^(-rate * maturity))
/// * `sqrt_maturity_sigma` - Square root of (maturity * volatility^2)
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let discount = 0.99;
/// let sigma = 0.3;
/// let maturity: f64 = 2.0;
/// let sqrt_maturity_sigma = sigma * maturity.sqrt();
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
/// # Arguments
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `sigma` - Volatility of the underlying asset
/// * `maturity` - Time to maturity in years
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma = 0.3;
/// let maturity = 1.0;
/// let call = black_scholes::call(stock, strike, rate, sigma, maturity);
/// ```
pub fn call(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    call_discount(s, k, (-rate * maturity).exp(), maturity.sqrt() * sigma)
}

/// Returns delta of a BS call option
///
/// Delta measures the rate of change of the option price with respect to changes in the underlying asset price.
///
/// # Arguments
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `sigma` - Volatility of the underlying asset
/// * `maturity` - Time to maturity in years
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma = 0.3;
/// let maturity = 1.0;
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
    } else if s > k {
        1.0
    } else {
        0.0
    }
}

/// Returns gamma of a BS call option
///
/// Gamma measures the rate of change in the delta with respect to changes in the underlying price.
///
/// # Arguments
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `sigma` - Volatility of the underlying asset
/// * `maturity` - Time to maturity in years
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma = 0.3;
/// let maturity = 1.0;
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
/// Vega measures sensitivity to volatility. Vega is the derivative of the option value with respect to the volatility of the underlying asset.
///
/// # Arguments
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `sigma` - Volatility of the underlying asset
/// * `maturity` - Time to maturity in years
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma = 0.3;
/// let maturity = 1.0;
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
/// Theta measures the sensitivity of the option price to time decay. It represents the rate of change of the option value with respect to time.
///
/// # Arguments
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `sigma` - Volatility of the underlying asset
/// * `maturity` - Time to maturity in years
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma = 0.3;
/// let maturity = 1.0;
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
/// Rho measures sensitivity to the interest rate. It is the derivative of the option value with respect to the risk-free interest rate.
///
/// # Arguments
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `sigma` - Volatility of the underlying asset
/// * `maturity` - Time to maturity in years
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma = 0.3;
/// let maturity = 1.0;
/// let rho = black_scholes::call_rho(
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
/// # Arguments
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `discount` - Discount factor (e^(-rate * maturity))
/// * `sqrt_maturity_sigma` - Square root of (maturity * volatility^2)
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let discount = 0.99;
/// let sigma = 0.3;
/// let maturity: f64 = 2.0;
/// let sqrt_maturity_sigma = sigma * maturity.sqrt();
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
/// # Arguments
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `sigma` - Volatility of the underlying asset
/// * `maturity` - Time to maturity in years
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma = 0.3;
/// let maturity = 1.0;
/// let put = black_scholes::put(stock, strike, rate, sigma, maturity);
/// ```
pub fn put(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    put_discount(s, k, (-rate * maturity).exp(), maturity.sqrt() * sigma)
}

/// Returns delta of a BS put option
///
/// Delta measures the rate of change of the option price with respect to changes in the underlying asset price.
///
/// # Arguments
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `sigma` - Volatility of the underlying asset
/// * `maturity` - Time to maturity in years
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma = 0.3;
/// let maturity = 1.0;
/// let delta = black_scholes::put_delta(
///     stock, strike, rate, sigma, maturity
/// );
/// ```
pub fn put_delta(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    let sqrt_maturity_sigma = maturity.sqrt() * sigma;
    if sqrt_maturity_sigma > 0.0 {
        let discount = (-rate * maturity).exp();
        let d1 = d1(s, k, discount, sqrt_maturity_sigma);
        cum_norm(d1) - 1.0
    } else if k > s {
        -1.0
    } else {
        0.0
    }
}

/// Returns gamma of a BS put option
///
/// Gamma measures the rate of change in the delta with respect to changes in the underlying price.
/// For put options, gamma is identical to call options.
///
/// # Arguments
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `sigma` - Volatility of the underlying asset
/// * `maturity` - Time to maturity in years
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma = 0.3;
/// let maturity = 1.0;
/// let gamma = black_scholes::put_gamma(
///     stock, strike, rate, sigma, maturity
/// );
/// ```
pub fn put_gamma(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    call_gamma(s, k, rate, sigma, maturity) //same as call
}

/// Returns vega of a BS put option
///
/// Vega measures sensitivity to volatility. For put options, vega is identical to call options.
///
/// # Arguments
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `sigma` - Volatility of the underlying asset
/// * `maturity` - Time to maturity in years
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma = 0.3;
/// let maturity = 1.0;
/// let vega = black_scholes::put_vega(
///     stock, strike, rate, sigma, maturity
/// );
/// ```
pub fn put_vega(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    call_vega(s, k, rate, sigma, maturity) //same as call
}

/// Returns theta of a BS put option
///
/// Theta measures the sensitivity of the option price to time decay. It represents the rate of change of the option value with respect to time.
///
/// # Arguments
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `sigma` - Volatility of the underlying asset
/// * `maturity` - Time to maturity in years
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma = 0.3;
/// let maturity = 1.0;
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
/// Rho measures sensitivity to the interest rate. It is the derivative of the option value with respect to the risk-free interest rate.
///
/// # Arguments
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `sigma` - Volatility of the underlying asset
/// * `maturity` - Time to maturity in years
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma = 0.3;
/// let maturity = 1.0;
/// let rho = black_scholes::put_rho(
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

/// Returns vanna of a BS call option
///
/// Vanna measures the sensitivity of delta to changes in volatility. It is the second-order derivative of the option price with respect to the underlying price and volatility.
///
/// # Arguments
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `sigma` - Volatility of the underlying asset
/// * `maturity` - Time to maturity in years
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma = 0.3;
/// let maturity = 1.0;
/// let vanna = black_scholes::call_vanna(
///     stock, strike, rate, sigma, maturity
/// );
/// ```
pub fn call_vanna(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    let sqrt_t = maturity.sqrt();
    let sqrt_maturity_sigma = sqrt_t * sigma;
    if sqrt_maturity_sigma > 0.0 {
        let discount = (-rate * maturity).exp();
        let d1 = d1(s, k, discount, sqrt_maturity_sigma);
        -inc_norm(d1) * (d1 - sqrt_maturity_sigma) / sigma
    } else {
        0.0
    }
}
/// Returns vanna of a BS put option
///
/// Vanna measures the sensitivity of delta to changes in volatility. For put options, vanna is identical to call options.
///
/// # Arguments
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `sigma` - Volatility of the underlying asset
/// * `maturity` - Time to maturity in years
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma = 0.3;
/// let maturity = 1.0;
/// let vanna = black_scholes::put_vanna(
///     stock, strike, rate, sigma, maturity
/// );
/// ```
pub fn put_vanna(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    call_vanna(s, k, rate, sigma, maturity)
}
/// Returns vomma of a BS call option
///
/// Vomma (Volga) measures the second-order sensitivity of the option price to changes in volatility. It is the second derivative of the option price with respect to volatility.
///
/// # Arguments
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `sigma` - Volatility of the underlying asset
/// * `maturity` - Time to maturity in years
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma = 0.3;
/// let maturity = 1.0;
/// let vomma = black_scholes::call_vomma(
///     stock, strike, rate, sigma, maturity
/// );
/// ```
pub fn call_vomma(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    let sqrt_t = maturity.sqrt();
    let sqrt_maturity_sigma = sqrt_t * sigma;
    if sqrt_maturity_sigma > 0.0 {
        let discount = (-rate * maturity).exp();
        let d1 = d1(s, k, discount, sqrt_maturity_sigma);
        let d2 = d1 - sqrt_maturity_sigma;
        s * inc_norm(d1) * d1 * d2 * maturity / (sqrt_maturity_sigma)
    } else {
        0.0
    }
}

/// Returns vomma of a BS put option
///
/// Vomma (Volga) measures the second-order sensitivity of the option price to changes in volatility. For put options, vomma is identical to call options.
///
/// # Arguments
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `sigma` - Volatility of the underlying asset
/// * `maturity` - Time to maturity in years
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma = 0.3;
/// let maturity = 1.0;
/// let vomma = black_scholes::put_vomma(
///     stock, strike, rate, sigma, maturity
/// );
/// ```
pub fn put_vomma(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    call_vomma(s, k, rate, sigma, maturity)
}

/// Returns charm of a BS call option
///
/// Charm (Delta decay) measures the rate of change of delta over time. It is the second-order derivative of the option price with respect to time and underlying price.
///
/// # Arguments
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `sigma` - Volatility of the underlying asset
/// * `maturity` - Time to maturity in years
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma = 0.3;
/// let maturity = 1.0;
/// let charm = black_scholes::call_charm(
///     stock, strike, rate, sigma, maturity
/// );
/// ```
pub fn call_charm(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    let sqrt_t = maturity.sqrt();
    let sqrt_maturity_sigma = sqrt_t * sigma;
    if sqrt_maturity_sigma > 0.0 {
        let discount = (-rate * maturity).exp();
        let d1 = d1(s, k, discount, sqrt_maturity_sigma);
        let d2 = d1 - sqrt_maturity_sigma;
        -inc_norm(d1) * (2.0 * rate * maturity - d2 * sqrt_maturity_sigma)
            / (2.0 * maturity * sqrt_maturity_sigma)
    } else {
        0.0 // TODO: check that this is true....
    }
}

/// Returns charm of a BS put option
///
/// Charm (Delta decay) measures the rate of change of delta over time. For put options, charm is identical to call options.
///
/// # Arguments
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `sigma` - Volatility of the underlying asset
/// * `maturity` - Time to maturity in years
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let rate = 0.05;
/// let sigma = 0.3;
/// let maturity = 1.0;
/// let charm = black_scholes::put_charm(
///     stock, strike, rate, sigma, maturity
/// );
/// ```
pub fn put_charm(s: f64, k: f64, rate: f64, sigma: f64, maturity: f64) -> f64 {
    call_charm(s, k, rate, sigma, maturity)
}

const SQRT_TWO_PI: f64 = 2.0 * SQRT_2 / FRAC_2_SQRT_PI;
//Corrado and Miller (1996)
fn approximate_vol(price: f64, s: f64, k: f64, rate: f64, maturity: f64) -> f64 {
    let discount = (-rate * maturity).exp();
    let x = k * discount;
    let coef = SQRT_TWO_PI / (s + x);
    let helper_1 = s - x;
    let c1 = price - helper_1 * 0.5;
    let c2 = c1.powi(2);
    let c3 = helper_1.powi(2) * FRAC_1_PI;
    let bridge_1 = c2 - c3;
    let bridge_m = bridge_1.max(0.0).sqrt();
    coef * (c1 + bridge_m) / maturity.sqrt()
}
/// Returns implied volatility from a call option with initial guess
///
/// Calculates the implied volatility that would result in the given option price using Newton-Raphson method.
///
/// # Arguments
/// * `price` - Market price of the option
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `maturity` - Time to maturity in years
/// * `initial_guess` - Initial guess for the volatility
///
/// # Returns
/// * `Ok(volatility)` if the algorithm converges
/// * `Err(error_code)` if the algorithm fails to converge
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
/// Calculates the implied volatility that would result in the given option price using Newton-Raphson method.
/// Uses an approximation formula for the initial guess.
///
/// # Arguments
/// * `price` - Market price of the option
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `maturity` - Time to maturity in years
///
/// # Returns
/// * `Ok(volatility)` if the algorithm converges
/// * `Err(error_code)` if the algorithm fails to converge
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
/// Calculates the implied volatility that would result in the given option price using Newton-Raphson method.
///
/// # Arguments
/// * `price` - Market price of the option
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `maturity` - Time to maturity in years
/// * `initial_guess` - Initial guess for the volatility
///
/// # Returns
/// * `Ok(volatility)` if the algorithm converges
/// * `Err(error_code)` if the algorithm fails to converge
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
/// Calculates the implied volatility that would result in the given option price using Newton-Raphson method.
/// Uses put-call parity to convert the put price to an equivalent call price for the initial guess.
///
/// # Arguments
/// * `price` - Market price of the option
/// * `s` - Current stock price
/// * `k` - Strike price
/// * `rate` - Risk-free interest rate
/// * `maturity` - Time to maturity in years
///
/// # Returns
/// * `Ok(volatility)` if the algorithm converges
/// * `Err(error_code)` if the algorithm fails to converge
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

/// Container for option prices and Greeks.
///
/// Contains all the calculated prices and Greeks for both call and put options.
#[derive(Debug, Serialize)]
pub struct PricesAndGreeks {
    /// Price of the call option
    pub call_price: f64,
    /// Delta of the call option
    pub call_delta: f64,
    /// Gamma of the call option
    pub call_gamma: f64,
    /// Theta of the call option
    pub call_theta: f64,
    /// Vega of the call option
    pub call_vega: f64,
    /// Rho of the call option
    pub call_rho: f64,
    /// Vanna of the call option
    pub call_vanna: f64,
    /// Vomma of the call option
    pub call_vomma: f64,
    /// Charm of the call option
    pub call_charm: f64,
    /// Price of the put option
    pub put_price: f64,
    /// Delta of the put option
    pub put_delta: f64,
    /// Gamma of the put option
    pub put_gamma: f64,
    /// Theta of the put option
    pub put_theta: f64,
    /// Vega of the put option
    pub put_vega: f64,
    /// Rho of the put option
    pub put_rho: f64,
    /// Vanna of the put option
    pub put_vanna: f64,
    /// Vomma of the put option
    pub put_vomma: f64,
    /// Charm of the put option
    pub put_charm: f64,
}
/// Returns call and put prices and Greeks.
/// Due to caching the complex computations
/// (such as N(d1)), this implementation is
/// faster if you need to obtain all the
/// information for a given stock price
/// and strike price.
///
/// # Arguments
/// * `stock` - Current stock price
/// * `strike` - Strike price
/// * `rate` - Risk-free interest rate
/// * `sigma` - Volatility of the underlying asset
/// * `maturity` - Time to maturity in years
///
/// # Returns
/// A `PricesAndGreeks` struct containing all prices and Greeks
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
// Not using annualised dividend yield (i.e. using Black-Scholes, not Merton formula):
// See https://github.com/danielhstahl/black_scholes_rust/issues/25 referring to https://www.macroption.com/black-scholes-formula/
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
        let call_vega = stock * pdf_d1 * sqrt_maturity;
        let call_rho = k_discount * maturity * cdf_d2;
        let call_vanna = call_vega / stock * (1.0 - d1 / sqrt_maturity_sigma);
        let call_vomma = call_vega * d1 * d2 / sigma;
        let call_charm = -pdf_d1 * (2.0 * rate * maturity - d2 * sqrt_maturity_sigma)
            / (2.0 * maturity * sqrt_maturity_sigma);
        let put_price = call_price + k_discount - stock;
        let put_delta = cdf_d1 - 1.0;
        let put_gamma = call_gamma;
        let put_theta =
            -stock * pdf_d1 * sigma / (2.0 * sqrt_maturity) + rate * k_discount * (1.0 - cdf_d2);
        let put_vega = call_vega;
        let put_rho = -1.0 * k_discount * maturity * (1.0 - cdf_d2);
        let put_vanna = call_vanna;
        let put_vomma = call_vomma;
        let put_charm = call_charm;
        PricesAndGreeks {
            call_price,
            call_delta,
            call_gamma,
            call_theta,
            call_vega,
            call_rho,
            call_vanna,
            call_vomma,
            call_charm,
            put_price,
            put_delta,
            put_gamma,
            put_theta,
            put_vega,
            put_rho,
            put_vanna,
            put_vomma,
            put_charm,
        }
    } else {
        PricesAndGreeks {
            call_price: max_or_zero(stock - strike),
            call_delta: if stock > strike { 1.0 } else { 0.0 },
            call_gamma: 0.0,
            call_theta: 0.0,
            call_vega: 0.0,
            call_rho: 0.0,
            call_vanna: 0.0,
            call_vomma: 0.0,
            call_charm: 0.0,
            put_price: max_or_zero(strike - stock),
            put_delta: if strike > stock { -1.0 } else { 0.0 },
            put_gamma: 0.0,
            put_theta: 0.0,
            put_vega: 0.0,
            put_rho: 0.0,
            put_vanna: 0.0,
            put_vomma: 0.0,
            put_charm: 0.0,
        }
    }
}

/// Returns call and put prices and Greeks using Black-Scholes-Merton formula.
///
/// The Black-Scholes-Merton model extends the Black-Scholes model to account for continuous dividends.
/// If `dividend_yield` is 0, this gives same results as `compute_all` using Black-Scholes formula
/// but `compute_all` will be slightly less compute intensive.
///
/// # Arguments
/// * `stock` - Current stock price
/// * `strike` - Strike price
/// * `sigma` - Volatility of the underlying asset
/// * `risk_free_rate` - Annualized continuously compounded risk-free interest rate
/// * `dividend_yield` - Annualized continuously compounded dividend yield
/// * `maturity` - Time to maturity in years
///
/// # Returns
/// A `PricesAndGreeks` struct containing all prices and Greeks
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let sigma = 0.3;
/// let risk_free_rate = 0.05;
/// let dividend_yield = 0.02;
/// let maturity = 1.0;
/// let all_prices_and_greeks = black_scholes::bsm_compute_all(
///     stock,
///     strike,
///     sigma,
///     risk_free_rate,
///     dividend_yield,
///     maturity,
/// );
/// ```
pub fn bsm_compute_all(
    stock: f64,
    strike: f64,
    sigma: f64,
    risk_free_rate: f64,
    dividend_yield: f64,
    maturity: f64,
) -> PricesAndGreeks {
    let dividend = (-dividend_yield * maturity).exp();
    let discount = (-risk_free_rate * maturity).exp();
    let sqrt_maturity = maturity.sqrt();
    let sqrt_maturity_sigma = sqrt_maturity * sigma;
    let k_discount = strike * discount;
    if sqrt_maturity_sigma > 0.0 {
        let d1 = ((stock / strike).ln()
            + (risk_free_rate - dividend_yield + 0.5 * sigma.powi(2)) * maturity)
            / sqrt_maturity_sigma;
        let d2 = d1 - sqrt_maturity_sigma;
        let cdf_d1 = cum_norm(d1);
        let cdf_d2 = cum_norm(d2);
        let pdf_d1 = inc_norm(d1);

        let call_price = stock * dividend * cdf_d1 - k_discount * cdf_d2;
        let call_delta = dividend * cdf_d1;
        let call_gamma = dividend * pdf_d1 / (stock * sqrt_maturity_sigma);
        let call_theta = -dividend * stock * pdf_d1 * sigma / (2.0 * sqrt_maturity)
            - risk_free_rate * k_discount * cdf_d2
            + dividend_yield * stock * dividend * cdf_d1;
        let call_vega = stock * pdf_d1 * sqrt_maturity;
        let call_rho = k_discount * maturity * cdf_d2;
        let call_vanna = call_vega / stock * (1.0 - d1 / sqrt_maturity_sigma);
        let call_vomma = call_vega * d1 * d2 / sigma;
        let charm_part = dividend
            * pdf_d1
            * (2.0 * (risk_free_rate - dividend_yield) * maturity - d2 * sqrt_maturity_sigma)
            / (2.0 * maturity * sqrt_maturity_sigma);
        let call_charm = dividend_yield * dividend * cdf_d1 - charm_part;

        let put_price = call_price + k_discount - stock * dividend;
        let put_delta = dividend * (cdf_d1 - 1.0);
        let put_gamma = call_gamma;
        let put_theta = -dividend * stock * pdf_d1 * sigma / (2.0 * sqrt_maturity)
            + risk_free_rate * k_discount * (1.0 - cdf_d2)
            - dividend_yield * stock * dividend * (1.0 - cdf_d1);
        let put_vega = call_vega;
        let put_rho = -k_discount * maturity * (1.0 - cdf_d2);
        let put_vanna = call_vanna;
        let put_vomma = call_vomma;
        let put_charm = -dividend_yield * dividend * (1.0 - cdf_d1) - charm_part;
        PricesAndGreeks {
            call_price,
            call_delta,
            call_gamma,
            call_theta,
            call_vega,
            call_rho,
            call_vanna,
            call_vomma,
            call_charm,
            put_price,
            put_delta,
            put_gamma,
            put_theta,
            put_vega,
            put_rho,
            put_vanna,
            put_vomma,
            put_charm,
        }
    } else {
        PricesAndGreeks {
            call_price: max_or_zero(stock - strike),
            call_delta: if stock > strike { 1.0 } else { 0.0 },
            call_gamma: 0.0,
            call_theta: 0.0,
            call_vega: 0.0,
            call_rho: 0.0,
            call_vanna: 0.0,
            call_vomma: 0.0,
            call_charm: 0.0,
            put_price: max_or_zero(strike - stock),
            put_delta: if strike > stock { -1.0 } else { 0.0 },
            put_gamma: 0.0,
            put_theta: 0.0,
            put_vega: 0.0,
            put_rho: 0.0,
            put_vanna: 0.0,
            put_vomma: 0.0,
            put_charm: 0.0,
        }
    }
}

/// Validates that the input parameters for option pricing are within reasonable bounds.
///
/// # Arguments
/// * `s` - Stock price (must be positive)
/// * `k` - Strike price (must be positive)
/// * `rate` - Risk-free rate (no bounds)
/// * `sigma` - Volatility (must be non-negative)
/// * `maturity` - Time to maturity (must be non-negative)
///
/// # Returns
/// `true` if all parameters are valid, `false` otherwise
pub fn validate_params(s: f64, k: f64, _rate: f64, sigma: f64, maturity: f64) -> bool {
    s > 0.0 && k > 0.0 && sigma >= 0.0 && maturity >= 0.0 && s.is_finite() && k.is_finite() && sigma.is_finite() && maturity.is_finite()
}

/// Returns call and put prices and Greeks using the Black 76 model for options on futures.
///
/// The Black 76 model is used for pricing options on futures contracts. It is similar to Black-Scholes
/// but uses the forward price instead of the spot price.
///
/// # Arguments
/// * `forward_price` - Forward price of the underlying asset
/// * `strike` - Strike price
/// * `rate` - Risk-free interest rate
/// * `sigma` - Volatility of the underlying asset
/// * `maturity` - Time to maturity in years
///
/// # Returns
/// A `PricesAndGreeks` struct containing all prices and Greeks
///
/// # Examples
///
/// ```
/// let forward_price = 55.0;
/// let strike = 50.0;
/// let rate = 0.05;
/// let sigma = 0.2;
/// let maturity = 1.0;
/// let all_prices_and_greeks = black_scholes::black76(
///     forward_price,
///     strike,
///     rate,
///     sigma,
///     maturity,
/// );
/// ```
pub fn black76(
    forward_price: f64,
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
        let ln_f_s = (forward_price / strike).ln();

        let d1 = (ln_f_s + 0.5 * sigma.powi(2) * maturity) / sqrt_maturity_sigma;
        let d2 = d1 - sqrt_maturity_sigma;
        let cdf_d1 = cum_norm(d1); // often noted `N(d1)`
        let cdf_d2 = cum_norm(d2);
        let pdf_d1 = inc_norm(d1); // often noted `n(d1)`

        let call_price = discount * (forward_price * cdf_d1 - strike * cdf_d2);

        let call_delta = cdf_d1 * discount;
        let call_gamma = discount * pdf_d1 / (forward_price * sqrt_maturity_sigma);
        let call_theta = -forward_price * discount * pdf_d1 * sigma / (2.0 * sqrt_maturity)
            - rate * k_discount * cdf_d2
            + rate * forward_price * discount * cdf_d1;
        let call_vega = forward_price * discount * pdf_d1 * sqrt_maturity;
        let call_rho = -maturity * discount * (forward_price * cdf_d1 - strike * cdf_d2);
        let call_vanna = (call_vega / forward_price) * (1.0 - d1 / sqrt_maturity_sigma);
        let call_vomma = call_vega * d1 * d2 / sigma;

        let charm_part = pdf_d1
            * ((sigma / (4.0 * sqrt_maturity)) - (ln_f_s / (2.0 * sqrt_maturity_sigma * maturity)));
        let call_charm = discount * ((-rate * cdf_d1) + charm_part);

        // Deduce Put price from Call price using the put-call parity: https://en.wikipedia.org/wiki/Put%E2%80%93call_parity
        //  `Call - Put = Discount . (Fwd - K)`
        //
        // Can also find put price from call price formula using `cdf(x) + cdf(-x) == 1` and `pdf(x) == pfd(-x)`
        let put_price = call_price + discount * (strike - forward_price);

        let put_delta = discount * (cdf_d1 - 1.0);
        let put_gamma = call_gamma;

        let put_theta = -forward_price * discount * pdf_d1 * sigma / (2.0 * sqrt_maturity)
            + rate * k_discount * (1.0 - cdf_d2)
            - rate * forward_price * discount * (1.0 - cdf_d1);
        let put_vega = call_vega;
        let put_rho =
            -maturity * discount * (strike * (1.0 - cdf_d2) - forward_price * (1.0 - cdf_d1));
        let put_vanna = call_vanna;
        let put_vomma = call_vomma;
        let put_charm = discount * ((rate * (1.0 - cdf_d1)) + charm_part);

        PricesAndGreeks {
            call_price,
            call_delta,
            call_gamma,
            call_theta,
            call_vega,
            call_rho,
            call_vanna,
            call_vomma,
            call_charm,
            put_price,
            put_delta,
            put_gamma,
            put_theta,
            put_vega,
            put_rho,
            put_vanna,
            put_vomma,
            put_charm,
        }
    } else {
        PricesAndGreeks {
            call_price: max_or_zero(forward_price - strike),
            call_delta: if forward_price > strike { 1.0 } else { 0.0 },
            call_gamma: 0.0,
            call_theta: 0.0,
            call_vega: 0.0,
            call_rho: 0.0,
            call_vanna: 0.0,
            call_vomma: 0.0,
            put_price: max_or_zero(strike - forward_price),
            put_delta: if strike > forward_price { -1.0 } else { 0.0 },
            put_gamma: 0.0,
            put_theta: 0.0,
            put_vega: 0.0,
            put_rho: 0.0,
            put_vanna: 0.0,
            put_vomma: 0.0,
            call_charm: 0.0,
            put_charm: 0.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::*;
    use rand::SeedableRng;
    use rand::distr::{Distribution, Uniform};
    use rand::rngs::StdRng;
    use std::f64::consts::PI;

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
        assert_abs_diff_eq!(SQRT_TWO_PI, (2.0 * PI).sqrt(), epsilon = 0.000000001);
    }
    #[test]
    fn constants_are_correct() {
        assert_approx_eq!(FRAC_1_SQRT_2PI, (2.0 * PI).sqrt().recip());
    }
    #[test]
    fn cum_norm_opposite() {
        fn check(x: f64) {
            assert_abs_diff_eq!(cum_norm(x) + cum_norm(-x), 1.0, epsilon = 0.000000001);
        }
        check(0.0);
        check(0.1);
        check(0.2);
        check(0.9);
        check(1.0);
        check(2.0);
        check(10.0);
    }
    #[test]
    fn inc_norm_opposite() {
        fn check(x: f64) {
            assert_abs_diff_eq!(inc_norm(x), inc_norm(-x), epsilon = 0.000000001);
        }
        check(0.0);
        check(0.1);
        check(0.2);
        check(0.9);
        check(1.0);
        check(2.0);
        check(10.0);
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
        assert_abs_diff_eq!(sigma, approx_vol, epsilon = 0.01);
    }
    #[test]
    fn call_iv_works_with_broad_set_of_numbers() {
        let seed: [u8; 32] = [2; 32];
        let mut rng_seed = get_rng_seed(seed);
        let uniform = Uniform::new(0.0f64, 1.0).unwrap();
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
        let _iv = call_iv(price, s, k, rate, maturity).unwrap();
    }
    #[test]
    fn put_iv_works() {
        let sigma = 0.2;
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
    fn put_iv_returns_err_if_no_possible_solution() {
        // Use a clearly invalid case - put price higher than max possible value
        let price = 1000.0;  // Very high price that's impossible
        let s = 100.0;       // Underlying price
        let k = 100.0;       // Strike price
        let rate = 0.05;     // Interest rate
        let maturity = 1.0;  // Time to maturity

        assert!(put_iv(price, s, k, rate, maturity).is_err());
    }

    #[test]
    fn call_discount_with_zero_sqrt_maturity_sigma() {
        let s = 5.0;
        let k = 4.5;
        let discount = 0.99;
        let sqrt_maturity_sigma = 0.0;
        let price = call_discount(s, k, discount, sqrt_maturity_sigma);
        assert_eq!(price, max_or_zero(s - k));
    }

    #[test]
    fn put_discount_with_zero_sqrt_maturity_sigma() {
        let s = 5.0;
        let k = 4.5;
        let discount = 0.99;
        let sqrt_maturity_sigma = 0.0;
        let price = put_discount(s, k, discount, sqrt_maturity_sigma);
        assert_eq!(price, max_or_zero(k - s));
    }

    #[test]
    fn call_gamma_with_zero_volatility() {
        let s = 5.0;
        let k = 4.5;
        let rate = 0.05;
        let sigma = 0.0;
        let maturity = 1.0;
        let gamma = call_gamma(s, k, rate, sigma, maturity);
        assert_eq!(gamma, 0.0);
    }

    #[test]
    fn put_gamma_with_zero_volatility() {
        let s = 5.0;
        let k = 4.5;
        let rate = 0.05;
        let sigma = 0.0;
        let maturity = 1.0;
        let gamma = put_gamma(s, k, rate, sigma, maturity);
        assert_eq!(gamma, 0.0);
    }

    #[test]
    fn call_vega_with_zero_volatility() {
        let s = 5.0;
        let k = 4.5;
        let rate = 0.05;
        let sigma = 0.0;
        let maturity = 1.0;
        let vega = call_vega(s, k, rate, sigma, maturity);
        assert_eq!(vega, 0.0);
    }

    #[test]
    fn put_vega_with_zero_volatility() {
        let s = 5.0;
        let k = 4.5;
        let rate = 0.05;
        let sigma = 0.0;
        let maturity = 1.0;
        let vega = put_vega(s, k, rate, sigma, maturity);
        assert_eq!(vega, 0.0);
    }

    #[test]
    fn call_theta_with_zero_volatility() {
        let s = 5.0;
        let k = 4.5;
        let rate = 0.05;
        let sigma = 0.0;
        let maturity = 1.0;
        let theta = call_theta(s, k, rate, sigma, maturity);
        assert_eq!(theta, 0.0);
    }

    #[test]
    fn put_theta_with_zero_volatility() {
        let s = 5.0;
        let k = 4.5;
        let rate = 0.05;
        let sigma = 0.0;
        let maturity = 1.0;
        let theta = put_theta(s, k, rate, sigma, maturity);
        assert_eq!(theta, 0.0);
    }

    #[test]
    fn call_rho_with_zero_volatility() {
        let s = 5.0;
        let k = 4.5;
        let rate = 0.05;
        let sigma = 0.0;
        let maturity = 1.0;
        let rho = call_rho(s, k, rate, sigma, maturity);
        assert_eq!(rho, 0.0);
    }

    #[test]
    fn put_rho_with_zero_volatility() {
        let s = 5.0;
        let k = 4.5;
        let rate = 0.05;
        let sigma = 0.0;
        let maturity = 1.0;
        let rho = put_rho(s, k, rate, sigma, maturity);
        assert_eq!(rho, 0.0);
    }

    #[test]
    fn call_vanna_with_zero_volatility() {
        let s = 5.0;
        let k = 4.5;
        let rate = 0.05;
        let sigma = 0.0;
        let maturity = 1.0;
        let vanna = call_vanna(s, k, rate, sigma, maturity);
        assert_eq!(vanna, 0.0);
    }

    #[test]
    fn put_vanna_with_zero_volatility() {
        let s = 5.0;
        let k = 4.5;
        let rate = 0.05;
        let sigma = 0.0;
        let maturity = 1.0;
        let vanna = put_vanna(s, k, rate, sigma, maturity);
        assert_eq!(vanna, 0.0);
    }

    #[test]
    fn call_vomma_with_zero_volatility() {
        let s = 5.0;
        let k = 4.5;
        let rate = 0.05;
        let sigma = 0.0;
        let maturity = 1.0;
        let vomma = call_vomma(s, k, rate, sigma, maturity);
        assert_eq!(vomma, 0.0);
    }

    #[test]
    fn put_vomma_with_zero_volatility() {
        let s = 5.0;
        let k = 4.5;
        let rate = 0.05;
        let sigma = 0.0;
        let maturity = 1.0;
        let vomma = put_vomma(s, k, rate, sigma, maturity);
        assert_eq!(vomma, 0.0);
    }

    #[test]
    fn call_charm_with_zero_volatility() {
        let s = 5.0;
        let k = 4.5;
        let rate = 0.05;
        let sigma = 0.0;
        let maturity = 1.0;
        let charm = call_charm(s, k, rate, sigma, maturity);
        assert_eq!(charm, 0.0);
    }

    #[test]
    fn put_charm_with_zero_volatility() {
        let s = 5.0;
        let k = 4.5;
        let rate = 0.05;
        let sigma = 0.0;
        let maturity = 1.0;
        let charm = put_charm(s, k, rate, sigma, maturity);
        assert_eq!(charm, 0.0);
    }

    #[test]
    fn charm_edge_case_check() {
        // Test the TODO comment in call_charm function
        let s = 550.88;
        let sigma = 0.37;
        let k = 510.0;
        let rate = 0.0;
        let maturity = -0.09; // negative maturity
        let charm = call_charm(s, k, rate, sigma, maturity);
        assert_eq!(charm, 0.0); // This verifies the TODO comment is correct
    }

    #[test]
    fn black76_with_zero_volatility() {
        let forward_price = 55.0;
        let k = 50.0;
        let maturity = 1.0;
        let sigma = 0.0;
        let rate = 0.0025;
        let result = black76(forward_price, k, rate, sigma, maturity);

        // With zero volatility, call price should be max(0, forward - strike)
        assert_eq!(result.call_price, (forward_price - k).max(0.0));
        // With zero volatility, put price should be max(0, strike - forward)
        assert_eq!(result.put_price, (k - forward_price).max(0.0));
    }

    #[test]
    fn black76_put_call_parity() {
        let forward_price = 55.0;
        let k = 50.0;
        let maturity = 1.0;
        let sigma = 0.15;
        let rate = 0.0025;
        let result = black76(forward_price, k, rate, sigma, maturity);

        // Put-call parity for futures options: call - put = discount * (forward - strike)
        let discount = (-rate * maturity).exp();
        let expected_parity = discount * (forward_price - k);
        let actual_parity = result.call_price - result.put_price;
        assert_abs_diff_eq!(actual_parity, expected_parity, epsilon = 1e-10);
    }

    #[test]
    fn compute_all_with_zero_volatility() {
        let s = 5.0;
        let k = 4.5;
        let rate = 0.05;
        let sigma = 0.0;
        let maturity = 1.0;
        let result = compute_all(s, k, rate, sigma, maturity);

        assert_eq!(result.call_price, (s - k).max(0.0));
        assert_eq!(result.put_price, (k - s).max(0.0));
        assert_eq!(result.call_delta, if s > k { 1.0 } else { 0.0 });
        assert_eq!(result.put_delta, if k > s { -1.0 } else { 0.0 });
        assert_eq!(result.call_gamma, 0.0);
        assert_eq!(result.put_gamma, 0.0);
        assert_eq!(result.call_theta, 0.0);
        assert_eq!(result.put_theta, 0.0);
        assert_eq!(result.call_vega, 0.0);
        assert_eq!(result.put_vega, 0.0);
        assert_eq!(result.call_rho, 0.0);
        assert_eq!(result.put_rho, 0.0);
        assert_eq!(result.call_vanna, 0.0);
        assert_eq!(result.put_vanna, 0.0);
        assert_eq!(result.call_vomma, 0.0);
        assert_eq!(result.put_vomma, 0.0);
        assert_eq!(result.call_charm, 0.0);
        assert_eq!(result.put_charm, 0.0);
    }

    #[test]
    fn bsm_compute_all_with_zero_dividend() {
        let s = 5.0;
        let k = 4.5;
        let sigma = 0.3;
        let rate = 0.05;
        let q = 0.0;  // Zero dividend yield
        let maturity = 1.0;

        let bs_result = compute_all(s, k, rate, sigma, maturity);
        let bsm_result = bsm_compute_all(s, k, sigma, rate, q, maturity);

        // Results should be nearly identical when dividend yield is zero
        assert_abs_diff_eq!(bs_result.call_price, bsm_result.call_price, epsilon = 1e-10);
        assert_abs_diff_eq!(bs_result.put_price, bsm_result.put_price, epsilon = 1e-10);
        assert_abs_diff_eq!(bs_result.call_delta, bsm_result.call_delta, epsilon = 1e-10);
        assert_abs_diff_eq!(bs_result.put_delta, bsm_result.put_delta, epsilon = 1e-10);
        assert_abs_diff_eq!(bs_result.call_gamma, bsm_result.call_gamma, epsilon = 1e-10);
        assert_abs_diff_eq!(bs_result.put_gamma, bsm_result.put_gamma, epsilon = 1e-10);
        assert_abs_diff_eq!(bs_result.call_theta, bsm_result.call_theta, epsilon = 1e-10);
        assert_abs_diff_eq!(bs_result.put_theta, bsm_result.put_theta, epsilon = 1e-10);
        assert_abs_diff_eq!(bs_result.call_vega, bsm_result.call_vega, epsilon = 1e-10);
        assert_abs_diff_eq!(bs_result.put_vega, bsm_result.put_vega, epsilon = 1e-10);
        assert_abs_diff_eq!(bs_result.call_rho, bsm_result.call_rho, epsilon = 1e-10);
        assert_abs_diff_eq!(bs_result.put_rho, bsm_result.put_rho, epsilon = 1e-10);
        assert_abs_diff_eq!(bs_result.call_vanna, bsm_result.call_vanna, epsilon = 1e-10);
        assert_abs_diff_eq!(bs_result.put_vanna, bsm_result.put_vanna, epsilon = 1e-10);
        assert_abs_diff_eq!(bs_result.call_vomma, bsm_result.call_vomma, epsilon = 1e-10);
        assert_abs_diff_eq!(bs_result.put_vomma, bsm_result.put_vomma, epsilon = 1e-10);
        assert_abs_diff_eq!(bs_result.call_charm, bsm_result.call_charm, epsilon = 1e-10);
        assert_abs_diff_eq!(bs_result.put_charm, bsm_result.put_charm, epsilon = 1e-10);
    }

    #[test]
    fn bsm_compute_all_with_positive_dividend() {
        let s = 100.0;
        let k = 100.0;
        let sigma = 0.2;
        let rate = 0.05;
        let q = 0.02;  // Positive dividend yield
        let maturity = 1.0;

        let bs_result = compute_all(s, k, rate, sigma, maturity);
        let bsm_result = bsm_compute_all(s, k, sigma, rate, q, maturity);

        // With positive dividend yield, call prices should be lower and put prices higher
        assert!(bsm_result.call_price < bs_result.call_price);
        assert!(bsm_result.put_price > bs_result.put_price);
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
    fn put_greeks_work_2() {
        let s = 150.0;
        let sigma = 0.37;
        let k = 160.0;
        let rate = 0.03;
        let maturity = 0.08;
        assert_abs_diff_eq!(
            put_delta(s, k, rate, sigma, maturity),
            -0.71,
            epsilon = 0.01
        );
        assert_abs_diff_eq!(
            put_theta(s, k, rate, sigma, maturity),
            -30.26,
            epsilon = 0.01
        );
        assert_abs_diff_eq!(put_vega(s, k, rate, sigma, maturity), 14.62, epsilon = 0.01);
        assert_abs_diff_eq!(put_rho(s, k, rate, sigma, maturity), -9.46, epsilon = 0.01);
        assert_abs_diff_eq!(
            put_gamma(s, k, rate, sigma, maturity),
            0.022,
            epsilon = 0.01
        );
        assert_abs_diff_eq!(
            put_vanna(s, k, rate, sigma, maturity),
            0.602,
            epsilon = 0.01
        );
        assert_abs_diff_eq!(
            put_charm(s, k, rate, sigma, maturity),
            -1.490,
            epsilon = 0.01
        );
        assert_abs_diff_eq!(
            put_vomma(s, k, rate, sigma, maturity),
            13.821,
            epsilon = 0.01
        );
    }

    #[test]
    fn call_greeks_works_2() {
        let s = 150.0;
        let sigma = 0.37;
        let k = 160.0;
        let rate = 0.03;
        let maturity = 0.08;
        assert_abs_diff_eq!(
            call_delta(s, k, rate, sigma, maturity),
            0.29,
            epsilon = 0.01
        );
        assert_abs_diff_eq!(
            call_theta(s, k, rate, sigma, maturity),
            -35.04,
            epsilon = 0.01
        );
        assert_abs_diff_eq!(
            call_vega(s, k, rate, sigma, maturity),
            14.62,
            epsilon = 0.01
        );
        assert_abs_diff_eq!(call_rho(s, k, rate, sigma, maturity), 3.31, epsilon = 0.01);
        assert_abs_diff_eq!(
            call_gamma(s, k, rate, sigma, maturity),
            0.022,
            epsilon = 0.01
        );
        assert_abs_diff_eq!(
            call_vanna(s, k, rate, sigma, maturity),
            0.602,
            epsilon = 0.01
        );
        assert_abs_diff_eq!(
            call_charm(s, k, rate, sigma, maturity),
            -1.490,
            epsilon = 0.01
        );
        assert_abs_diff_eq!(
            call_vomma(s, k, rate, sigma, maturity),
            13.821,
            epsilon = 0.01
        );
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
            call_vanna: call_vanna_result,
            call_vomma: call_vomma_result,
            call_charm: call_charm_result,
            put_price,
            put_delta,
            put_gamma,
            put_theta,
            put_vega,
            put_rho,
            ..
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

        assert_approx_eq!(call_vanna_result, call_vanna(s, k, rate, sigma, maturity));
        assert_approx_eq!(call_vomma_result, call_vomma(s, k, rate, sigma, maturity));
        assert_approx_eq!(call_charm_result, call_charm(s, k, rate, sigma, maturity));
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

    #[test]
    fn bsm_compute_all_zero_q() {
        // If `q` is zero, `bsm_compute_all`` is equiv. to `compute_all`
        let s = 550.88;
        let sigma = 0.37;
        let k = 510.0;
        let rate = 0.05;
        let q = 0.0;
        let maturity = 0.09;

        let r0 = compute_all(s, k, rate, sigma, maturity);
        let r1 = bsm_compute_all(s, k, sigma, rate, q, maturity);

        let PricesAndGreeks {
            call_price,
            call_delta,
            call_gamma,
            call_theta,
            call_vega,
            call_rho,
            call_vanna,
            call_vomma,
            call_charm,
            put_price,
            put_delta,
            put_gamma,
            put_theta,
            put_vega,
            put_rho,
            put_vanna,
            put_vomma,
            put_charm,
        } = r0;
        macro_rules! check {
            ($field:ident) => {{
                assert_approx_eq!($field, r1.$field);
            }};
        }
        check!(call_price);
        check!(call_delta);
        check!(call_gamma);
        check!(call_theta);
        check!(call_vega);
        check!(call_rho);
        check!(call_vanna);
        check!(call_vomma);
        check!(call_charm);
        check!(put_price);
        check!(put_delta);
        check!(put_gamma);
        check!(put_theta);
        check!(put_vega);
        check!(put_rho);
        check!(put_vanna);
        check!(put_vomma);
        check!(put_charm);
    }
    #[test]
    fn bsm_compute_all_works() {
        // Compare value to https://github.com/CarloLepelaars/blackscholes
        // (with fix for theta: https://github.com/CarloLepelaars/blackscholes/pull/26)
        let s = 550.88;
        let sigma = 0.37;
        let k = 510.0;
        let rate = 0.05;
        let q = 0.03;
        let maturity = 0.09;

        let PricesAndGreeks {
            call_price,
            call_delta,
            call_gamma,
            call_theta,
            call_vega,
            call_rho,
            call_vanna,
            call_vomma,
            call_charm,
            put_price,
            put_delta,
            put_gamma,
            put_theta,
            put_vega,
            put_rho,
            put_vanna,
            put_vomma,
            put_charm,
        } = bsm_compute_all(s, k, sigma, rate, q, maturity);
        assert_approx_eq!(call_price, 49.9003280);
        assert_approx_eq!(call_delta, 0.7761726197638565);
        assert_approx_eq!(call_gamma, 0.004850905458223078);
        assert_approx_eq!(call_theta, -106.82167402360267);
        assert_approx_eq!(call_vega, 49.15340973000121);
        assert_approx_eq!(call_rho, 33.99098802331468);
        assert_approx_eq!(call_vanna, -0.5268153918879713);
        assert_approx_eq!(call_vomma, 66.72271336536006);
        assert_approx_eq!(call_charm, 1.049818266353188);
        assert_approx_eq!(put_price, 8.215853933440584);
        assert_approx_eq!(put_delta, -0.22113102195785656);
        assert_approx_eq!(put_gamma, 0.004850905458223078);
        assert_approx_eq!(put_theta, -97.91800512749833);
        assert_approx_eq!(put_vega, 49.15340973000121);
        assert_approx_eq!(put_rho, -11.702926017862612);
        assert_approx_eq!(put_vanna, -0.5268153918879713);
        assert_approx_eq!(put_vomma, 66.72271336536006);
        assert_approx_eq!(put_charm, 1.0198991571015366);
    }
    #[test]
    fn black76_works() {
        let s = 55.;
        let k = 50.0;
        let maturity = 1.0;
        let sigma = 0.15;
        let rate = 0.0025;
        let PricesAndGreeks {
            call_price,
            call_delta,
            call_gamma,
            call_theta,
            call_vega,
            call_rho,
            call_vanna,
            call_vomma,
            call_charm,
            put_price,
            put_delta,
            put_gamma,
            put_theta,
            put_vega,
            put_rho,
            put_vanna,
            put_vomma,
            put_charm,
        } = black76(s, k, rate, sigma, maturity);
        assert_approx_eq!(call_price, 6.234516);
        assert_approx_eq!(call_delta, 0.759371);
        assert_approx_eq!(call_gamma, 0.037478);
        assert_approx_eq!(call_theta, -1.259855);
        assert_approx_eq!(call_vega, 17.005889);
        assert_approx_eq!(call_rho, -6.234516);
        assert_approx_eq!(call_vanna, -1.155166);
        assert_approx_eq!(call_vomma, 45.134728);
        assert_approx_eq!(call_charm, -0.088535); // value not verified externally :(

        assert_approx_eq!(put_price, 1.247001);
        assert_approx_eq!(put_delta, -0.238131);
        assert_approx_eq!(put_gamma, 0.037478);
        assert_approx_eq!(put_theta, -1.272324);
        assert_approx_eq!(put_vega, 17.005889);
        assert_approx_eq!(put_rho, -1.247001);
        assert_approx_eq!(put_vanna, -1.155166);
        assert_approx_eq!(put_vomma, 45.134728);
        assert_approx_eq!(put_charm, -0.086042); // value not verified externally :(
    }

    #[test]
    fn call_delta_with_negative_maturity_works() {
        let s = 550.88;
        let sigma = 0.37;
        let k = 510.0;
        let rate = 0.0;
        let maturity = -0.09;
        assert_approx_eq!(call_delta(s, k, rate, sigma, maturity), 1.0);
    }

    #[test]
    fn put_delta_with_negative_maturity_works() {
        let s = 550.88;
        let sigma = 0.37;
        let k = 510.0;
        let rate = 0.0;
        let maturity = -0.09;
        assert_approx_eq!(put_delta(s, k, rate, sigma, maturity), 0.0);
    }

    #[test]
    fn call_delta_with_negative_maturity_and_s_less_than_k_works() {
        let s = 510.0;
        let sigma = 0.37;
        let k = 550.88;
        let rate = 0.0;
        let maturity = -0.09;
        assert_approx_eq!(call_delta(s, k, rate, sigma, maturity), 0.0);
    }

    #[test]
    fn put_delta_with_negative_maturity_and_s_less_than_k_works() {
        let s = 510.0;
        let sigma = 0.37;
        let k = 550.88;
        let rate = 0.0;
        let maturity = -0.09;
        assert_approx_eq!(put_delta(s, k, rate, sigma, maturity), -1.0);
    }

    #[test]
    fn call_gamma_with_negative_maturity_works() {
        let s = 550.88;
        let sigma = 0.37;
        let k = 510.0;
        let rate = 0.0;
        let maturity = -0.09;
        assert_approx_eq!(call_gamma(s, k, rate, sigma, maturity), 0.0);
    }

    #[test]
    fn theta_with_negative_maturity_works() {
        let s = 550.88;
        let sigma = 0.37;
        let k = 510.0;
        let rate = 0.0;
        let maturity = -0.09;
        assert_approx_eq!(call_theta(s, k, rate, sigma, maturity), 0.0);
        assert_approx_eq!(put_theta(s, k, rate, sigma, maturity), 0.0);
    }

    #[test]
    fn rho_with_negative_maturity_works() {
        let s = 550.88;
        let sigma = 0.37;
        let k = 510.0;
        let rate = 0.0;
        let maturity = -0.09;
        assert_approx_eq!(call_rho(s, k, rate, sigma, maturity), 0.0);
        assert_approx_eq!(put_rho(s, k, rate, sigma, maturity), 0.0);
    }

    #[test]
    fn vanna_with_negative_maturity_works() {
        let s = 550.88;
        let sigma = 0.37;
        let k = 510.0;
        let rate = 0.0;
        let maturity = -0.09;
        assert_approx_eq!(call_vanna(s, k, rate, sigma, maturity), 0.0);
        assert_approx_eq!(put_vanna(s, k, rate, sigma, maturity), 0.0);
    }

    #[test]
    fn vomma_with_negative_maturity_works() {
        let s = 550.88;
        let sigma = 0.37;
        let k = 510.0;
        let rate = 0.0;
        let maturity = -0.09;
        assert_approx_eq!(call_vomma(s, k, rate, sigma, maturity), 0.0);
        assert_approx_eq!(put_vomma(s, k, rate, sigma, maturity), 0.0);
    }

    #[test]
    fn charm_with_negative_maturity_works() {
        let s = 550.88;
        let sigma = 0.37;
        let k = 510.0;
        let rate = 0.0;
        let maturity = -0.09;
        assert_approx_eq!(call_charm(s, k, rate, sigma, maturity), 0.0); // TODO: check that this is true....
        assert_approx_eq!(put_charm(s, k, rate, sigma, maturity), 0.0);
    }
}
