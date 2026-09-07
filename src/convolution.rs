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
use crate::algebra::{Monoid, Ring, Semiring};

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

/// The monoid algebra `R[C]` truncated to `[0, n)`, as a semiring with values `Vec<R::Value>`.
///
/// # Definition
/// `R[C]` is the free `R`-module on the monoid `C` with the product extended bilinearly from `C`,
/// i.e. the convolution. Its zero is the zero vector, its one is `δ_id`, addition is pointwise, and
/// multiplication is `C::convolve`.
///
/// # Contract
/// The truncation to `[0, n)` is a ring homomorphism (see [`Convolution`]).
///
/// # Complexity
/// - Space: O(n) per value
pub struct MonoidAlgebra<R, C> {
    ring: R,
    conv: C,
    len: usize,
}

impl<R: Semiring, C: Convolution<R>> MonoidAlgebra<R, C> {
    /// # Panics
    /// Panics if `id >= len`.
    pub fn new(ring: R, conv: C, len: usize) -> Self {
        assert!(
            conv.id() < len,
            "id must be less than len: id={}, len={len}",
            conv.id()
        );
        Self { ring, conv, len }
    }
    pub fn len(&self) -> usize {
        self.len
    }
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}
impl<R: Semiring, C: Convolution<R>> Semiring for MonoidAlgebra<R, C> {
    type Value = Vec<R::Value>;
    fn zero(&self) -> Vec<R::Value> {
        (0..self.len).map(|_| self.ring.zero()).collect()
    }
    fn one(&self) -> Vec<R::Value> {
        let mut e = self.zero();
        e[self.conv.id()] = self.ring.one();
        e
    }
    /// # Panics
    /// Panics if the lengths of `a`, `b` and `self` differ.
    fn add(&self, a: &Vec<R::Value>, b: &Vec<R::Value>) -> Vec<R::Value> {
        assert!(
            a.len() == self.len && b.len() == self.len,
            "length mismatch: lhs={}, rhs={}, len={}",
            a.len(),
            b.len(),
            self.len
        );
        a.iter().zip(b).map(|(x, y)| self.ring.add(x, y)).collect()
    }
    /// # Panics
    /// Panics if the lengths of `a`, `b` and `self` differ.
    fn mul(&self, a: &Vec<R::Value>, b: &Vec<R::Value>) -> Vec<R::Value> {
        assert!(
            a.len() == self.len && b.len() == self.len,
            "length mismatch: lhs={}, rhs={}, len={}",
            a.len(),
            b.len(),
            self.len
        );
        self.conv.convolve(&self.ring, a, b)
    }
}
impl<R: Ring, C: Convolution<R>> Ring for MonoidAlgebra<R, C> {
    /// # Panics
    /// Panics if the length of `a` differs from `self`.
    fn neg(&self, a: &Vec<R::Value>) -> Vec<R::Value> {
        assert!(
            a.len() == self.len,
            "length mismatch: a={}, len={}",
            a.len(),
            self.len
        );
        a.iter().map(|x| self.ring.neg(x)).collect()
    }
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
