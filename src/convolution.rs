pub mod and;
pub mod xor;

use crate::algebra::{Monoid, Semiring};

/// A transform diagonalizing `R[Self]`, where `Self` is a monoid structure on `usize`.
///
/// # Contract
/// For `f`, `g` of the same supported length,
/// - (additivity) `T(f + g) = T(f) + T(g)`
/// - (multiplicativity) `T(f * g) = T(f) .* T(g)`, where `*` is the product of `R[Self]`
///   and `.*` is the pointwise product.
pub trait Transform<R: Semiring>: Monoid<Value = usize> {
    /// # Panics
    /// Panics if `f.len()` is not a supported length.
    fn transform(&self, ring: &R, f: &mut [R::Value]);
}
/// # Contract
/// `inverse . transform = id` on every supported length.
pub trait InverseTransform<R: Semiring>: Transform<R> {
    /// # Panics
    /// Panics if `f.len()` is not a supported length.
    fn inverse_transform(&self, ring: &R, f: &mut [R::Value]);
}

/// Returns `f * g` in `R[Self]`, computed as `T^{-1}(T(f) .* T(g))`.
///
/// # Complexity
/// - Time: 2 transform + 1 inverse transform + O(n)
/// - Space: 2 transform + 1 inverse transform + O(1)
///
/// # Panics
/// Panics if `f.len() != g.len()`.
/// Panics if the length is not supported by `t`.
pub fn convolve<R: Semiring, T: Transform<R> + InverseTransform<R>>(
    ring: &R,
    t: &T,
    mut f: Vec<R::Value>,
    mut g: Vec<R::Value>,
) -> Vec<R::Value> {
    assert_eq!(
        f.len(),
        g.len(),
        "length must agree: {} != {}",
        f.len(),
        g.len()
    );
    t.transform(ring, &mut f);
    t.transform(ring, &mut g);
    for (a, b) in f.iter_mut().zip(g.iter()) {
        *a = ring.mul(a, b);
    }
    t.inverse_transform(ring, &mut f);
    f
}
