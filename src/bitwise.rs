//! The bitwise operations `and`, `or` and `xor`, the set functions on `[n]`, and the algebras they
//! form.
//!
//! A subset `S` of `[n] = {0, ..., n - 1}` is identified with the integer in `[0, 2^n)` whose bit
//! `v` is set iff `v` is in `S`, so that `and`, `or` and `xor` are the intersection, the union and
//! the symmetric difference. A set function `f` assigns a value of `R` to each subset of `[n]`, and
//! of that length, and `S`, `T`, `U` denote subsets of `[n]`.

pub mod and;
pub mod or;
pub mod set_power_series;
pub mod xor;
