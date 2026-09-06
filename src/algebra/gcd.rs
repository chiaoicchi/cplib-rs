use crate::algebra::{Monoid, Zero};
use crate::num::gcd::gcd;

/// The gcd monoid of the non-negative integers.
///
/// # Definition
/// `gcd` is associative, commutative and idempotent, so `(N, 0, gcd)` is a meet-semilattice under
/// the divisibility order `x <= y` iff `x | y`; its top is `0`, since every `x` divides `0`, and
/// the top is the identity of the meet: `gcd(0, x) = x`.
///
/// Under `x -> (v_p(x))_p`, the positive integers with `gcd` are the direct product of the chains
/// `(N, min)` over the primes `p`, with finitely many nonzero coordinates.
///
/// # Contract
/// `%` is the Euclidian remainder on `T`, and all values are non-negative. For signed integers the
/// sign of `gcd(a, b)` depends on the signs of the inputs, and the monoid laws hold only up to
/// sign.
pub struct Gcd<T>(std::marker::PhantomData<T>);
impl<T> Gcd<T> {
    pub const fn new() -> Self {
        Self(std::marker::PhantomData)
    }
}
impl<T> Default for Gcd<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Clone for Gcd<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for Gcd<T> {}

impl<T: Clone + PartialEq + std::ops::Rem<Output = T> + Zero> Monoid for Gcd<T> {
    type Value = T;
    fn id(&self) -> T {
        T::zero()
    }
    fn op(&self, a: &T, b: &T) -> T {
        gcd(a.clone(), b.clone())
    }
}
