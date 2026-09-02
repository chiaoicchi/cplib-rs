pub mod additive;
pub mod closures;

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

/// A group trait
///
/// # Definition
/// A 4-tuple `(Value, id, op, inv)` is called a group if it satisfies:
/// - (monoid) `(Value, id, op)` forms a monoid.
/// - (inverse) `op(a, inv(a)) == op(inv(a), a) == id()` for all `a`.
pub trait Group: Monoid {
    fn inv(&self, a: &Self::Value) -> Self::Value;
}

/// An action trait
///
/// # Definition
/// A map `act: U x T -> T` is called an external law of composition on `T`,
/// where the elements of `U` are called operators.
pub trait Action<T, U> {
    fn act(&self, f: &U, x: &T) -> T;
}

/// A marker trait for commutative operations.
///
/// # Definition
/// An operation `op` is called commutative if `op(a, b) == op(b, a)` for all `a`, `b`.
pub trait Commutative {}
