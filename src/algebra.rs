pub mod additive;
pub mod canonical;
pub mod closures;

/// A monoid trait
///
/// # Definition
/// A triple `(Value, id, op)` is called a monoid if it satisfies:
/// - (associativity) `op(op(a, b), c) = op(a, op(b, c))` for all `a`, `b`, `c`.
/// - (identity) `op(id(), a) = op(a, id()) = a` for all `a`.
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
/// - (inverse) `op(a, inv(a)) = op(inv(a), a) = id()` for all `a`.
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

/// A semiring trait
///
/// # Definition
/// A 5-tuple `(Value, zero, one, add, mul)` is called a semiring if it satisfies:
/// - (additive monoid) `(Value, zero, add)` forms a commutative monoid.
/// - (multiplicative monoid) `(Value, one, mul)` forms a monoid.
/// - (distributivity) `mul(add(a, b), c) = add(mul(a, c), mul(b, c))` and
///   `mul(a, add(b, c)) = add(mul(a, b), mul(a, c))` for all `a`, `b`, `c`.
/// - (annihilation) `mul(zero(), a) = mul(a, zero()) = zero()` for all `a`.
pub trait Semiring {
    type Value;
    fn zero(&self) -> Self::Value;
    fn one(&self) -> Self::Value;
    fn add(&self, a: &Self::Value, b: &Self::Value) -> Self::Value;
    fn mul(&self, a: &Self::Value, b: &Self::Value) -> Self::Value;
}

/// A ring trait
///
/// # Definition
/// A 6-tuple `(Value, zero, one, add, neg, mul)` is called a ring if it satisfies:
/// - (additive group) `(Value, zero, add, neg)` forms a commutative group.
/// - (multiplicative monoid) `(Value, one, mul)` forms a monoid.
/// - (distributivity) `mul(add(a, b), c) = add(mul(a, c), mul(b, c))` and
///   `mul(a, add(b, c)) = add(mul(a, b), mul(a, c))` for all `a`, `b`, `c`.
pub trait Ring: Semiring {
    fn neg(&self, a: &Self::Value) -> Self::Value;
}

/// A field trait
///
/// # Definition
/// A 7-tuple `(Value, zero, one, add, neg, mul, inv)` is called a field if it satisfies:
/// - (additive group) `(Value, zero, add, neg)` forms a commutative group.
/// - (multiplicative group) `(Value\{zero}, one, mul, inv)` forms a commutative group.
/// - (distributivity) `mul(add(a, b), c) = add(mul(a, c), mul(b, c))` and
///   `mul(a, add(b, c)) = add(mul(a, b), mul(a, c))` for all `a`, `b`, `c`.
///
/// # Contract
/// `inv(zero())` is not defined; implementations may panic.
pub trait Field: Ring {
    fn inv(&self, a: &Self::Value) -> Self::Value;
}

/// A marker trait for commutative operations.
///
/// # Definition
/// An operation `op` is called commutative if `op(a, b) = op(b, a)` for all `a`, `b`.
pub trait Commutative {}

/// A type with a distinguished element `zero`.
///
/// # Contract
/// `zero` is the additive identity of `T`, i.e. `T::zero() + a = a + T::zero() == a` for all `a`.
pub trait Zero {
    fn zero() -> Self;
}

/// A type with a distinguished element `one`.
///
/// # Contract
/// `one` is the multiplicative identity of `T`, i.e. `T::one() * a = a * T::one() = a` for all `a`.
pub trait One {
    fn one() -> Self;
}

/// A type with a distinguished function `inv`.
///
/// # Contract
/// `inv` returns the multiplicative inverse in `T`, i.e. `a * inv(a) = inv(a) * a = one` for all
/// `a`.
pub trait Inv {
    type Output;
    fn inv(&self) -> Self::Output;
}
