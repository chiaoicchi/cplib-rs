//! Formal power series over a semiring `R`, truncated at a precision, as slices of coefficients.
//!
//! A series `f` is given by a slice with `f(k)`, the coefficient of `x^k`, at index `k`, and
//! `f(k) = 0` for `k >= f.len()`. Functions take a precision `n` and return the `n` coefficients of
//! the result in `R[[x]]/(x^n)`, which depends only on `f(0), ..., f(n - 1)`. Throughout this
//! module `n` is the precision.

use crate::algebra::RootOfUnity;
use crate::poly::poly_convolve;

pub mod elementary;
pub mod relaxed;

/// The product of `f` and `g` in `R[[x]]/(x^n)`, as its `n` coefficients.
///
/// # Definition
/// `(fg)(k) = Σ_{i+j=k} f(i) g(j)` for `k` in `[0, n)`.
///
/// # Complexity
/// - Time: O(n log n)
/// - Space: O(n)
///
/// # Panics
/// Panics if `a = min(f.len(), n)` and `b = min(g.len(), n)` are both greater than 32 and `R` has
/// no primitive `N`-th root of unity, where `N` is the least power of two at least `a + b - 1`.
pub fn fps_convolve<R: RootOfUnity>(
    ring: &R,
    mut f: Vec<R::Value>,
    mut g: Vec<R::Value>,
    n: usize,
) -> Vec<R::Value> {
    f.truncate(n);
    g.truncate(n);
    let mut h = poly_convolve(ring, f, g);
    h.resize_with(n, || ring.zero());
    h
}
