| [Linux][lin-link] |  [Coveralls][cov-link]  |
| :---------------: | :-------------------: |
| ![lin-badge]      | ![cov-badge]          |

[lin-badge]: https://github.com/danielhstahl/black_scholes_rust/workflows/Rust/badge.svg
[lin-link]: https://github.com/danielhstahl/black_scholes_rust/actions

[cov-badge]: https://coveralls.io/repos/github/danielhstahl/black_scholes_rust/badge.svg?branch=master
[cov-link]:  https://coveralls.io/repos/github/danielhstahl/black_scholes_rust


# black_scholes_rust

This is a simple Black Scholes option calculator written in rust.  Documentation is on [docs.rs](https://docs.rs/black_scholes).

## breaking changes

The move from 0.4 to 0.5 results changed the IV api to return a `Result<f64, f64>` rather than an `f64`.

## using black_scholes_rust
Put the following in your Cargo.toml:

```toml
[dependencies]
black_scholes = "0.10.1"
```

Import and use:

```rust
extern crate black_scholes;
let stock = 5.0;
let strike = 4.5;
let rate = 0.01;
let discount = 0.99;
let sigma = 0.3;
let maturity = 2.0;
let sqrt_maturity_sigma = sigma*maturity.sqrt();
let price = black_scholes::call_discount(
    stock, strike, discount,
    sqrt_maturity_sigma
);
//or 
let price = black_scholes::call(
    stock, strike, rate,
    sigma, maturity
);
```

