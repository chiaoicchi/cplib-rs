use crate::algebra::{Bounded, Commutative, Idempotent, Monoid, Semigroup};

/// The min monoid of a totally ordered set.
///
/// # Definition
/// `min` is associative, commutative and idempotent, so `(T, min)` is a meet-semilattice under
/// `Ord`, and the greatest element of `T` is the identity of the meet: `min(max_value, x) = x`.
pub struct Min<T>(std::marker::PhantomData<T>);
impl<T> Min<T> {
    pub const fn new() -> Self {
        Self(std::marker::PhantomData)
    }
}
impl<T> Default for Min<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Clone for Min<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for Min<T> {}

impl<T: Clone + Ord> Semigroup for Min<T> {
    type Value = T;
    fn op(&self, a: &T, b: &T) -> T {
        T::min(a.clone(), b.clone())
    }
}
impl<T: Clone + Ord + Bounded> Monoid for Min<T> {
    fn id(&self) -> T {
        T::min_value()
    }
}
impl<T: Clone + Ord> Commutative for Min<T> {}
impl<T: Clone + Ord> Idempotent for Min<T> {}
