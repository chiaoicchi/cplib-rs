use crate::algebra::{Action, Group, Monoid, Semigroup};

/// A monoid built from closures.
///
/// # Contract
/// `(T, id, op)` must form a monoid.
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

/// A group built from closures.
///
/// # Contract
/// `(T, id, op, inv)` must form a group.
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

/// An action built from closures.
pub struct FnAction<F> {
    pub act: F,
}
impl<T, U, F: Fn(&U, &T) -> T> Action<T, U> for FnAction<F> {
    fn act(&self, f: &U, x: &T) -> T {
        (self.act)(f, x)
    }
}
