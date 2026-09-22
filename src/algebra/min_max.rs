use crate::algebra::{Bounded, Commutative, Idempotent, Monoid, Semigroup};

/// The monoid of `T` under `min`.
///
/// # Definition
/// `op(a, b) = min(a, b)`, and `id()` is `T::max_value()`.
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
        T::max_value()
    }
}
impl<T: Clone + Ord> Commutative for Min<T> {}
impl<T: Clone + Ord> Idempotent for Min<T> {}

/// The monoid of `T` under `max`.
///
/// # Definition
/// `op(a, b) = max(a, b)`, and `id()` is `T::min_value()`.
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
