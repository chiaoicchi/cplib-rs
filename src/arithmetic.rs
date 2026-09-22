//! Arithmetic functions on `N` and the algebras they form.
//!
//! An arithmetic function `f` assigns a value of `R` to each `x ∈ N = {0, 1, ...}`, and is stored
//! as a slice of length `n` with `f(x)` at index `x`. Throughout this module `n` is that length,
//! and `x`, `y`, `d` denote elements of `[0, n)`.

pub mod dirichlet;
pub mod gcd;
pub mod lcm;
pub mod lpf;
pub mod multiplicative;
pub mod prime;
