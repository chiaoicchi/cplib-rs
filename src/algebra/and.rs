use crate::algebra::{Commutative, Idempotent, Monoid, Semigroup, Zero};

/// The and monoid of `(F_2)^n`.
///
/// # Definition
/// `({0, 1}, 1, &)` is a commutative idempotent monoid, i.e. a meet-semilattice with top `1`.
/// Its `n`-fold direct product has the operation defined componentwise:
/// - `(x & y)_i = x_i & y_i`. Identifying `x` in `[0, 2^n)` with its binary digits
///   `(x_0, ..., x_{n-1})`, this is the bitwise and, with identity `2^n - 1`.
///
/// # Contract
/// `id()` is `!T::zero()`, the all-ones value of `T`, which lies outside `[0, 2^n)`; `[0, 2^n)` is
/// a subsemigroup whose own identity is `2^n - 1`. Convolution and transforms use only `op`, so
/// this does not affect them.
pub struct And<T>(std::marker::PhantomData<T>);
impl<T> And<T> {
    pub const fn new() -> Self {
        Self(std::marker::PhantomData)
    }
}
impl<T> Default for And<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Clone for And<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for And<T> {}

impl<T: Clone + std::ops::BitAnd<Output = T> + std::ops::Not<Output = T>> Semigroup for And<T> {
    type Value = T;
    fn op(&self, a: &T, b: &T) -> T {
        a.clone() & b.clone()
    }
}
impl<T: Clone + std::ops::BitAnd<Output = T> + std::ops::Not<Output = T> + Zero> Monoid for And<T> {
    fn id(&self) -> T {
        !T::zero()
    }
}
impl<T: Clone + std::ops::BitOr<Output = T>> Commutative for And<T> {}
impl<T: Clone + std::ops::BitOr<Output = T>> Idempotent for And<T> {}
