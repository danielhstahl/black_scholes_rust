| [Linux][lin-link] |  [Codecov][cov-link]  |
| :---------------: | :-------------------: |
| ![lin-badge]      | ![cov-badge]          |

[lin-badge]: https://travis-ci.org/phillyfan1138/black_scholes_rust.svg?branch=master "Travis build status"
[lin-link]:  https://travis-ci.org/phillyfan1138/black_scholes_rust "Travis build status"
[cov-badge]: https://codecov.io/gh/phillyfan1138/black_scholes_rust/branch/master/graph/badge.svg
[cov-link]:  https://codecov.io/gh/phillyfan1138/black_scholes_rust


# black_scholes_rust

This is a simple Black Scholes option calculator written in rust.  Documentation is on [crates.io](https://docs.rs/black_scholes/0.3.0/black_scholes/).

## using black_scholes_rust

```toml
[dependencies]
black_scholes = "0.3"
```

```rust
extern crate black_scholes;
let stock = 5.0;
let strike = 4.5;
let discount = 0.99;
let sigma = 0.3;
let maturity:f64 = 2.0;
let sqrt_maturity_sigma = sigma*maturity.sqrt();
let price = black_scholes::call_discount(
    stock, strike, discount, 
    sqrt_maturity_sigma
);
```

# tests

Note that the greeks don't have direct tests.  Please feel free to add tests.  However, these formulas are essentially and extensively tested in the [fang oost option](https://github.com/phillyfan1138/fang_oost_option_rust) library.  