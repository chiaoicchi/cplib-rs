//! Structures given by closures.

use crate::algebra::{Action, Group, Monoid, Semigroup};

/// A monoid given by closures.
///
/// # Definition
/// `op(a, b) = (self.op)(a, b)` and `id() = self.id`.
///
/// # Contract
/// `op` and `id` satisfy the laws of [`Monoid`].
pub struct FnMonoid<T, F> {
    pub id: T,
    pub op: F,
}
impl<T: Clone, F: Fn(&T, &T) -> T> Semigroup for FnMonoid<T, F> {
    type Value = T;
    fn op(&self, a: &T, b: &T) -> T {
        (self.op)(a, b)
    }
}
impl<T: Clone, F: Fn(&T, &T) -> T> Monoid for FnMonoid<T, F> {
    fn id(&self) -> T {
        self.id.clone()
    }
}

/// A group given by closures.
///
/// # Definition
/// `op(a, b) = (self.op)(a, b)`, `id() = self.id` and `inv(a) = (self.inv)(a)`.
///
/// # Contract
/// `op`, `id` and `inv` satisfy the laws of [`Group`].
pub struct FnGroup<T, F, G> {
    pub id: T,
    pub op: F,
    pub inv: G,
}
impl<T: Clone, F: Fn(&T, &T) -> T, G: Fn(&T) -> T> Semigroup for FnGroup<T, F, G> {
    type Value = T;
    fn op(&self, a: &T, b: &T) -> T {
        (self.op)(a, b)
    }
}
impl<T: Clone, F: Fn(&T, &T) -> T, G: Fn(&T) -> T> Monoid for FnGroup<T, F, G> {
    fn id(&self) -> T {
        self.id.clone()
    }
}
impl<T: Clone, F: Fn(&T, &T) -> T, G: Fn(&T) -> T> Group for FnGroup<T, F, G> {
    fn inv(&self, a: &T) -> T {
        (self.inv)(a)
    }
}

/// An action given by closures.
///
/// # Definition
/// `act(f, x) = (self.act)(f, x)`.
pub struct FnAction<F> {
    pub act: F,
}
impl<T, U, F: Fn(&U, &T) -> T> Action<T, U> for FnAction<F> {
    fn act(&self, f: &U, x: &T) -> T {
        (self.act)(f, x)
    }
}
