//! Exponential generating functions over a field `R`, truncated at a precision, as slices of
//! coefficients.
//!
//! A series `f` is given by a slice with `f(k)`, the coefficient of `x^k / k!`, at index `k`, and
//! `f(k) = 0` for `k >= f.len()`. Functions take a precision `n` and return the `n` coefficients of
//! the result. Throughout this module `n` is the precision.

use crate::algebra::RootOfUnity;
use crate::combinatorics::factorial::Factorial;
use crate::fps::fps_convolve;

/// The product of `f` and `g` as exponential generating functions, as its `n` coefficients.
///
/// # Definition
/// `(fg)(k) = Σ_{i+j=k} binomial(k, i) f(i) g(j)` for `k` in `[0, n)`.
///
/// # Complexity
/// - Time: O(n log n)
/// - Space: O(n)
///
/// # Panics
/// Panics if `a = min(f.len(), n)` and `b = min(g.len(), n)` are both greater than 32 and `R` has
/// no primitive `N`-th root of unity, where `N` is the least power of two at least `a + b - 1`.
pub fn egf_fps_convolve<R: RootOfUnity<Value: Clone> + Clone>(
    ring: &R,
    mut f: Vec<R::Value>,
    mut g: Vec<R::Value>,
    n: usize,
) -> Vec<R::Value> {
    let factorial = Factorial::new(ring.clone(), n);
    f.truncate(n);
    g.truncate(n);
    for (i, f) in f.iter_mut().enumerate() {
        *f = ring.mul(f, &factorial.inv_factorial(i));
    }
    for (i, g) in g.iter_mut().enumerate() {
        *g = ring.mul(g, &factorial.inv_factorial(i));
    }
    let mut h = fps_convolve(ring, f, g, n);
    for (i, h) in h.iter_mut().enumerate() {
        *h = ring.mul(h, &factorial.factorial(i));
    }
    h
}
