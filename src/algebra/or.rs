use crate::algebra::{Monoid, Zero};

/// The xor monoid of `(F_2)^n`.
///
/// # Definition
/// `F_2 = ({0, 1}, 0, |)` is a commutative monoid. `(F_2)^n` is the `n`-fold direct product of
/// `F_2`, whose operation is defined componentwise:
/// - `(x + y)_i = x_i + y_i` in `F_2`. Identifying `x` in `[0, 2^n)` with its binary digits
///   `(x_0, ..., x_{n-1})`.
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
