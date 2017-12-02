//! # black_scholes
//! A Black Scholes option pricing library.
extern crate special;
//pub use black_scholes::call;
////pub mod black_scholes {
fn cum_norm(x:f64)->f64 {
    //is this some weird scope issue?
    use special::Error;
    use std::f64::consts::SQRT_2;
    (x/SQRT_2).erf()*0.5+0.5
}
/// Returns BS call option formula.
///
/// # Examples
///
/// ```
/// let stock = 5.0;
/// let strike = 4.5;
/// let discount=(-0.05 as f64).exp();
/// let sigma=0.3;
/// assert_eq!(0.9848721043419868, black_scholes::call(stock, strike, discount, sigma));
/// ```
pub fn call(s:f64, k:f64, discount:f64, sigma:f64)->f64{
    if sigma>0.0{
        let d1=(s/(k*discount)).ln()/sigma+0.5*sigma;
        return s*cum_norm(d1)-k*discount*cum_norm(d1-sigma);
    }
    else{
        return if s>k {s-k} else {0.0};
    }
}
//}

#[cfg(test)]
mod tests {
    use super::*;
    //use super::black_scholes;
    #[test]
    fn call_formula_works() {
        assert_eq!(call(5.0, 4.5, (-0.05 as f64).exp(), 0.3), 0.9848721043419868);
    }
}//
//pub mod black_scholes;