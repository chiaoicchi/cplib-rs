//! Polynomials over a semiring `R`, as slices of coefficients.
//!
//! A polynomial `f` is stored with `f(k)`, the coefficient of `x^k`, at index `k`. Trailing zeros
//! are allowed, and the empty slice is the zero polynomial. Results are not normalized: their
//! length is determined by the lengths of the operands alone. Throughout this module `n` and `m`
//! are the lengths of the operands `f` and `g`.

pub mod calculus;
pub mod evaluation;
pub mod geometric;
pub mod interpolation;
pub mod subproduct_tree;
pub mod taylor_shift;

use crate::algebra::RootOfUnity;
use crate::periodic_function::cyclic_convolve;

/// The length of the shorter operand up to which [`poly_convolve`] multiplies directly.
const NAIVE_LIMIT: usize = 32;
/// The product of `f` and `g` in `R[x]`.
///
/// # Definition
/// `(fg)(k) = Σ_{i+j=k} f(i) g(j)`, the product of the monoid algebra `R[(N, +)] = R[x]`. The
/// result has length `n + m - 1`.
///
/// # Complexity
/// - Time: O((n + m) log (n + m))
/// - Space: O(n + m)
///
/// # Panics
/// Panics if both `f` and `g` are longer than 32 and `R` has no primitive `N`-th root of unity,
/// where `N` is the least power of two at least `n + m - 1`.
pub fn poly_convolve<R: RootOfUnity>(
    ring: &R,
    mut f: Vec<R::Value>,
    mut g: Vec<R::Value>,
) -> Vec<R::Value> {
    if f.is_empty() || g.is_empty() {
        return Vec::new();
    }
    let len = f.len() + g.len() - 1;
    if f.len().min(g.len()) <= NAIVE_LIMIT {
        let mut h: Vec<R::Value> = (0..len).map(|_| ring.zero()).collect();
        for (i, fi) in f.iter().enumerate() {
            for (hj, gj) in h[i..].iter_mut().zip(g.iter()) {
                *hj = ring.add(hj, &ring.mul(fi, gj));
            }
        }
        return h;
    }
    let n = len.next_power_of_two();
    f.resize_with(n, || ring.zero());
    g.resize_with(n, || ring.zero());
    let mut h = cyclic_convolve(ring, f, g);
    h.truncate(len);
    h
}
