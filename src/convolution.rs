pub mod and;
pub mod cyclic;
pub mod gcd;
pub mod lcm;
pub mod or;
pub mod xor;

use crate::algebra::cyclic::Cyclic;
use crate::algebra::{Monoid, RootOfUnity, Semiring};

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

/// Returns the product of `f` and `g` in `R[t]`, of length `f.len() + g.len() - 1`.
///
/// # Definition
/// `R[N] = R[t]` embeds into `R[Z/nZ] = R[t]/(t^n - 1)`, injectively on polynomials of degree less
/// than `n`. For `n >= f.len() + g.len() - 1` the product has degree less than `n`, so it is
/// recovered from the cyclic convolution of length `n` without wrap-around.
///
/// # Complexity
/// - Time: O(n log n) with `n` the least power of two at least `f.len() + g.len() - 1`
///
/// # Panics
/// Panics if `R` has no primitive `n`-th root of unity.
pub fn convolve_poly<R: RootOfUnity>(
    ring: &R,
    mut f: Vec<R::Value>,
    mut g: Vec<R::Value>,
) -> Vec<R::Value> {
    if f.is_empty() || g.is_empty() {
        return vec![];
    }
    let len = f.len() + g.len() - 1;
    let n = len.next_power_of_two();
    f.resize_with(n, || ring.zero());
    g.resize_with(n, || ring.zero());
    let mut h = convolve(ring, &Cyclic::new(n), f, g);
    h.truncate(len);
    h
}
