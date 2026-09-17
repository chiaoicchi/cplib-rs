//! Set functions on `[n]` and the algebras they form.
//!
//! A set function `f` assigns a value of `R` to each subset of `[n] = {0, ..., n - 1}`, and is
//! stored as a slice of length `2^n` with `f(S)` at the index whose bit `v` is set iff `v` is in
//! `S`. Throughout this module `n` is `log2` of that length, and `S`, `T`, `U` denote subsets of
//! `[n]`.

pub mod set_power_series;
pub mod subset;
pub mod superset;
pub mod xor;
