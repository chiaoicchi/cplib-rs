use crate::algebra::{Monoid, One, Zero};
use crate::num::gcd::lcm;

/// The lcm monoid of the non-negative integers.
///
/// # Definition
/// `lcm` is associative, commutative and idempotent, so `(N, 1, lcm)` is a join-semilattice under
/// the divisibility order `x <= y` iff `x | y`; its bottom is `1`, which is the identity of the
/// join: `lcm(1, x) = x`. Its top is `0`, which is absorbing: `lcm(0, x) = 0`.
///
/// Under `x -> (v_p(x))_p`, the positive integers with `lcm` are the direct product of the chains
/// `(N, max)` over the primes `p`, with finitely many nonzero coordinates.
///
/// # Contract
/// `%` is the Euclidean remainder on `T`, all values are non-negative, and every `lcm` fits in `T`.
pub struct Lcm<T>(std::marker::PhantomData<T>);
impl<T> Lcm<T> {
    pub const fn new() -> Self {
        Self(std::marker::PhantomData)
    }
}
impl<T> Default for Lcm<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Clone for Lcm<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for Lcm<T> {}

impl<
    T: Clone
        + PartialEq
        + std::ops::Mul<Output = T>
        + std::ops::Div<Output = T>
        + std::ops::Rem<Output = T>
        + Zero
        + One,
> Monoid for Lcm<T>
{
    type Value = T;
    fn id(&self) -> T {
        T::one()
    }
    fn op(&self, a: &T, b: &T) -> T {
        lcm(a.clone(), b.clone())
    }
}
