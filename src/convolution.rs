pub mod and;
pub mod cyclic;
pub mod dirichlet;
pub mod gcd;
pub mod lcm;
pub mod or;
pub mod xor;

use crate::algebra::and::And;
use crate::algebra::cyclic::Cyclic;
use crate::algebra::gcd::Gcd;
use crate::algebra::lcm::Lcm;
use crate::algebra::or::Or;
use crate::algebra::xor::Xor;
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

/// A monoid structure on `usize` whose monoid algebra `R[Self]` has a computable product.
///
/// # Contract
/// For `f`, `g` of the same length `n`, `convolve(ring, f, g)` is the product `f * g` in
/// `R[Self]` projected to `[0, n)`:
/// `(f * g)(z) = Σ_{op(x, y)=z} f(x)g(y)` over `x, y < n`, the terms with `op(x, y) >= n` being
/// dropped. The projection is a ring homomorphism when `{z >= n}` spans an ideal,
/// i.e. when `op(x, y) >= n` whenever `x >= n` or `y >= n`.
///
/// # Panics
/// Panics if `f.len() != g.len()`, or if the length is not supported by `Self`.
pub trait Convolution<R: Semiring>: Monoid<Value = usize> {
    fn convolve(&self, ring: &R, f: &[R::Value], g: &[R::Value]) -> Vec<R::Value>;
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
pub fn convolve_by_transform<R: Semiring, T: Transform<R> + InverseTransform<R>>(
    t: &T,
    ring: &R,
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

macro_rules! impl_convolution_by_transform {
    ($($t:ty),* $(,)?) => {$(
        impl<R: Semiring> Convolution<R> for $t
        where
            Self: Transform<R> + InverseTransform<R>,
            R::Value: Clone,
        {
            fn convolve(&self, ring: &R, f: &[R::Value], g: &[R::Value]) -> Vec<R::Value> {
                let (f, g) = (f.to_vec(), g.to_vec());
                convolve_by_transform(self, ring, f, g)
            }
        }
    )*};
}
impl_convolution_by_transform!(
    Xor<usize>,
    And<usize>,
    Or<usize>,
    Gcd<usize>,
    Lcm<usize>,
    Cyclic,
);
