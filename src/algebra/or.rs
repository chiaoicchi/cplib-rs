use crate::algebra::{Monoid, Zero};

/// The or monoid of `({0, 1}, 0, |)^n`.
///
/// # Definition
/// `({0, 1}, 0, |)` is a commutative idempotent monoid, i.e. a join-semilattice with bottom `0`.
/// Its `n`-fold direct product has the operation defined componentwise:
/// - `(x | y)_i = x_i | y_i`. Identifying `x` in `[0, 2^n)` with its binary digits
///   `(x_0, ..., x_{n-1})`, this is the bitwise or, with identity `0`.
pub struct Or<T>(std::marker::PhantomData<T>);
impl<T> Or<T> {
    pub const fn new() -> Self {
        Self(std::marker::PhantomData)
    }
}
impl<T> Default for Or<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Clone for Or<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for Or<T> {}

impl<T: Clone + std::ops::BitOr<Output = T> + Zero> Monoid for Or<T> {
    type Value = T;
    fn id(&self) -> T {
        T::zero()
    }
    fn op(&self, a: &T, b: &T) -> T {
        a.clone() | b.clone()
    }
}
