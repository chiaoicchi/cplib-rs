use crate::algebra::{Action, Group, Monoid, Ring, Semiring};

/// A monoid built from closures.
///
/// # Contract
/// `(T, id, op)` must form a monoid.
pub struct FnMonoid<T, F> {
    pub id: T,
    pub op: F,
}
impl<T: Clone, F: Fn(&T, &T) -> T> Monoid for FnMonoid<T, F> {
    type Value = T;
    fn id(&self) -> T {
        self.id.clone()
    }
    fn op(&self, a: &T, b: &T) -> T {
        (self.op)(a, b)
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
impl<T: Clone, F: Fn(&T, &T) -> T, G: Fn(&T) -> T> Monoid for FnGroup<T, F, G> {
    type Value = T;
    fn id(&self) -> T {
        self.id.clone()
    }
    fn op(&self, a: &T, b: &T) -> T {
        (self.op)(a, b)
    }
}
impl<T: Clone, F: Fn(&T, &T) -> T, G: Fn(&T) -> T> Group for FnGroup<T, F, G> {
    fn inv(&self, a: &T) -> T {
        (self.inv)(a)
    }
}

/// A semiring built from closures.
///
/// # Contract
/// `(T, zero, one, add, mul)` must form a semiring.
pub struct FnSemiring<T, F> {
    pub zero: T,
    pub one: T,
    pub add: F,
    pub mul: F,
}
impl<T: Clone, F: Fn(&T, &T) -> T> Semiring for FnSemiring<T, F> {
    type Value = T;
    fn zero(&self) -> T {
        self.zero.clone()
    }
    fn one(&self) -> T {
        self.one.clone()
    }
    fn add(&self, a: &T, b: &T) -> T {
        (self.add)(a, b)
    }
    fn mul(&self, a: &T, b: &T) -> T {
        (self.mul)(a, b)
    }
}

/// A ring built from closures.
///
/// # Contract
/// `()`
pub struct FnRing<T, F, G> {
    pub zero: T,
    pub one: T,
    pub add: F,
    pub mul: F,
    pub neg: G,
}
impl<T: Clone, F: Fn(&T, &T) -> T, G: Fn(&T) -> T> Semiring for FnRing<T, F, G> {
    type Value = T;
    fn zero(&self) -> T {
        self.zero.clone()
    }
    fn one(&self) -> T {
        self.one.clone()
    }
    fn add(&self, a: &T, b: &T) -> T {
        (self.add)(a, b)
    }
    fn mul(&self, a: &T, b: &T) -> T {
        (self.mul)(a, b)
    }
}
impl<T: Clone, F: Fn(&T, &T) -> T, G: Fn(&T) -> T> Ring for FnRing<T, F, G> {
    fn neg(&self, a: &T) -> T {
        (self.neg)(a)
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
