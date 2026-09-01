pub mod additive;

/// A monoid trait
///
/// # Definition
/// A triple `(Value, id, op)` is called a monoid if it satisfies:
/// - (associativity) `op(op(a, b), c) == op(a, op(b, c))` for all `a`, `b`, `c`.
/// - (identity) `op(id(), a) == op(a, id()) == a` for all `a`.
pub trait Monoid {
    type Value;
    fn id(&self) -> Self::Value;
    fn op(&self, a: &Self::Value, b: &Self::Value) -> Self::Value;
}

/// A monoid built from closures.
///
/// # Contract
/// `(T, id, op)` must form a monoid,
/// i.e. `op` is associative and `id` is its identity.
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

/// A group trait
///
/// # Definition
/// A 4-tuple `(Value, id, op, inv)` is called a group if it satisfies:
/// - (monoid) `(Value, id, op)` forms a monoid.
/// - (inverse) `op(a, inv(a)) == op(inv(a), a) == id()` for all `a`.
pub trait Group: Monoid {
    fn inv(&self, a: &Self::Value) -> Self::Value;
}

/// A group built from closures.
///
/// # Contract
/// `(T, id, op, inv)` must form a group,
/// i.e. `(T, id, op)` forms a monoid and `inv` maps each element to its inverse.
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

/// An action trait
///
/// # Definition
/// A map `act: U x T -> T` is called an external law of composition on `T`,
/// where the elements of `Map` are called operators.
pub trait Action<T, U> {
    fn act(&self, f: &U, x: &T) -> T;
}

/// An action built from closures.
///
/// # Contract
///
pub struct FnAction<F> {
    pub act: F,
}
impl<T, U, F: Fn(&U, &T) -> T> Action<T, U> for FnAction<F> {
    fn act(&self, f: &U, x: &T) -> T {
        (self.act)(f, x)
    }
}

/// A marker trait for commutative operations.
///
/// # Definition
/// An operation `op` is called commutative if `op(a, b) == op(b, a)` for all `a`, `b`.
pub trait Commutative {}
