/// A monoid trait
///
/// # Definition
/// ## Monoid
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
