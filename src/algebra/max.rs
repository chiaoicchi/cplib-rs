use crate::algebra::{Bounded, Commutative, Idempotent, Monoid, Semigroup};

/// The max monoid of a totally ordered set.
///
/// # Definition
/// `max` is associative, commutative and idempotent, so `(T, max)` is a join-semilattice under `Ord`,
/// and the least element of `T` is the identity of the join: `max(min_value, x) = x`.
pub struct Max<T>(std::marker::PhantomData<T>);
impl<T> Max<T> {
    pub const fn new() -> Self {
        Self(std::marker::PhantomData)
    }
}
impl<T> Default for Max<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Clone for Max<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for Max<T> {}

impl<T: Clone + Ord> Semigroup for Max<T> {
    type Value = T;
    fn op(&self, a: &T, b: &T) -> T {
        T::max(a.clone(), b.clone())
    }
}
impl<T: Clone + Ord + Bounded> Monoid for Max<T> {
    fn id(&self) -> T {
        T::min_value()
    }
}
impl<T: Clone + Ord> Commutative for Max<T> {}
impl<T: Clone + Ord> Idempotent for Max<T> {}
