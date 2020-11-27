| [Linux][lin-link] |  [Codecov][cov-link]  |
| :---------------: | :-------------------: |
| ![lin-badge]      | ![cov-badge]          |

[lin-badge]: https://github.com/danielhstahl/black_scholes_rust/workflows/Rust/badge.svg
[lin-link]: https://github.com/danielhstahl/black_scholes_rust/actions

[cov-badge]: https://codecov.io/gh/danielhstahl/black_scholes_rust/branch/master/graph/badge.svg
[cov-link]:  https://codecov.io/gh/danielhstahl/black_scholes_rust


# black_scholes_rust

This is a simple Black Scholes option calculator written in rust.  Documentation is on [docs.rs](https://docs.rs/black_scholes).

## breaking changes

The move from 0.4 to 0.5 results changed the IV api to return a `Result<f64, f64>` rather than an `f64`.

## using black_scholes_rust
Put the following in your Cargo.toml:

```toml
[dependencies]
black_scholes = "0.5"
```

Import and use:

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

Note that the greeks don't have direct tests.  Please feel free to add tests.  However, these formulas are essentially and extensively tested in the [fang oost option](https://github.com/danielhstahl/fang_oost_option_rust) library.