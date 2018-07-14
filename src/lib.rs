//! # black_scholes
//! A Black Scholes option pricing library.
extern crate special;
extern crate nrfind;
use std::f64::consts::PI;
use std::f64::consts::SQRT_2;
use special::Error;
#[macro_use]
#[cfg(test)]
extern crate approx;

fn cum_norm(x:f64)->f64 {
    (x/SQRT_2).erf()*0.5+0.5
}
fn inc_norm(x:f64)->f64 {
    (-x.powi(2)/2.0).exp()/(PI.sqrt()*SQRT_2)
}
/// Returns BS call option formula.
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
pub fn call(s:f64, k:f64, rate:f64, sigma:f64, maturity:f64)->f64{
    let sqrt_maturity_sigma=maturity.sqrt()*sigma;
    if sqrt_maturity_sigma>0.0{
        let discount=(-rate*maturity).exp();
        let d1=(s/(k*discount)).ln()/sqrt_maturity_sigma+0.5*sqrt_maturity_sigma;
        s*cum_norm(d1)-k*discount*cum_norm(d1-sqrt_maturity_sigma)
    }
    else{
        if s>k {s-k} else {0.0}
    }
}

pub fn call_delta(s:f64, k:f64, rate:f64, sigma:f64, maturity:f64)->f64{
    let sqrt_maturity_sigma=maturity.sqrt()*sigma;
    if sqrt_maturity_sigma>0.0{
        let discount=(-rate*maturity).exp();
        let d1=(s/(k*discount)).ln()/sqrt_maturity_sigma+0.5*sqrt_maturity_sigma;
        cum_norm(d1)
    }
    else{
        if s>k {1.0} else {0.0}
    }
}


pub fn call_gamma(s:f64, k:f64, rate:f64, sigma:f64, maturity:f64)->f64{
    let sqrt_maturity_sigma=maturity.sqrt()*sigma;    
    if sqrt_maturity_sigma>0.0{
        let discount=(-rate*maturity).exp();
        let d1=(s/(k*discount)).ln()/sqrt_maturity_sigma+0.5*sqrt_maturity_sigma;
        inc_norm(d1)/(s*sqrt_maturity_sigma)
    }
    else{
        0.0
    }
}

pub fn call_vega(s:f64, k:f64, rate:f64, sigma:f64, maturity:f64)->f64{
    let sqrt_maturity_sigma=maturity.sqrt()*sigma;    
    if sqrt_maturity_sigma>0.0{
        let discount=(-rate*maturity).exp();
        let d1=(s/(k*discount)).ln()/sqrt_maturity_sigma+0.5*sqrt_maturity_sigma;
        s*inc_norm(d1)*sqrt_maturity_sigma/sigma
    }
    else{
        0.0
    }
}
pub fn call_theta(s:f64, k:f64, rate:f64, sigma:f64, maturity:f64)->f64{
    let sqrt_t=maturity.sqrt();
    let sqrt_maturity_sigma=sqrt_t*sigma;    
    if sqrt_maturity_sigma>0.0{
        let discount=(-rate*maturity).exp();
        let d1=(s/(k*discount)).ln()/sqrt_maturity_sigma+0.5*sqrt_maturity_sigma;
        -s*inc_norm(d1)*sigma/(2.0*sqrt_t)-rate*k*discount*cum_norm(d1-sqrt_maturity_sigma)
    }
    else{
        0.0
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
pub fn put(s:f64, k:f64, rate:f64, sigma:f64, maturity:f64)->f64{
    let sqrt_maturity_sigma=maturity.sqrt()*sigma;
    if sqrt_maturity_sigma>0.0{
        let discount=(-rate*maturity).exp();  
        let d1=(s/(k*discount)).ln()/sqrt_maturity_sigma+0.5*sqrt_maturity_sigma;
        k*discount*cum_norm(sqrt_maturity_sigma-d1)-s*cum_norm(-d1)
    }
    else{
        if k>s {k-s} else {0.0}
    }
}

pub fn put_delta(s:f64, k:f64, rate:f64, sigma:f64, maturity:f64)->f64{
    let sqrt_maturity_sigma=maturity.sqrt()*sigma;
    if sqrt_maturity_sigma>0.0{
        let discount=(-rate*maturity).exp();  
        let d1=(s/(k*discount)).ln()/sqrt_maturity_sigma+0.5*sqrt_maturity_sigma;
        return cum_norm(d1)-1.0;
    }
    else{
        return if k>s {-1.0} else {0.0};
    }
}
pub fn put_gamma(s:f64, k:f64, rate:f64, sigma:f64, maturity:f64)->f64{
    call_gamma(s, k, rate, sigma, maturity)//same as call
}

pub fn put_vega(s:f64, k:f64, rate:f64, sigma:f64, maturity:f64)->f64{
    call_vega(s, k, rate, sigma, maturity) //same as call
}

pub fn put_theta(s:f64, k:f64, rate:f64, sigma:f64, maturity:f64)->f64{
    let sqrt_t=maturity.sqrt();
    let sqrt_maturity_sigma=sqrt_t*sigma;    
    if sqrt_maturity_sigma>0.0{
        let discount=(-rate*maturity).exp();
        let d1=(s/(k*discount)).ln()/sqrt_maturity_sigma+0.5*sqrt_maturity_sigma;
        -s*inc_norm(d1)*sigma/(2.0*sqrt_t)+rate*k*discount*cum_norm(-d1+sqrt_maturity_sigma)
    }
    else{
        0.0
    }
}

pub fn call_iv(price:f64, s:f64, k:f64, rate:f64, maturity:f64, initial_guess:f64)->f64{
    let obj_fn=|sigma|call(s, k, rate, sigma, maturity)-price;
    let dfn=|sigma|call_vega(s, k, rate, sigma, maturity);
    let precision=0.00000001;
    let iterations=20;
    nrfind::find_root(&obj_fn, &dfn, initial_guess, precision, iterations).unwrap()
}
pub fn put_iv(price:f64, s:f64, k:f64, rate:f64, maturity:f64, initial_guess:f64)->f64{
    let obj_fn=|sigma|put(s, k, rate, sigma, maturity)-price;
    let dfn=|sigma|put_vega(s, k, rate, sigma, maturity);
    let precision=0.00000001;
    let iterations=20;
    nrfind::find_root(&obj_fn, &dfn, initial_guess, precision, iterations).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
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
    fn call_iv_works(){
        let sigma=0.2;
        let initial_guess=0.5;
        let s=5.0;
        let k=4.5;
        let rate=0.05;
        let maturity=1.0;
        let price=call(s, k, rate, sigma, maturity);
        assert_abs_diff_eq!(call_iv(price, s, k, rate, maturity , initial_guess), sigma, epsilon=0.00000001);
    }
    #[test]
    fn put_iv_works(){
        let sigma=0.2;
        let initial_guess=0.5;
        let s=5.0;
        let k=4.5;
        let rate=0.05;
        let maturity=1.0;
        let price=put(s, k, rate, sigma, maturity);
        assert_abs_diff_eq!(put_iv(price, s, k, rate, maturity , initial_guess), sigma, epsilon=0.00000001);
    }
}
