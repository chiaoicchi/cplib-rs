pub mod canonical;
pub mod closures;
pub mod min_max;
pub mod monoid_algebra;
pub mod power;

/// A semigroup on `Value`.
///
/// # Contract
/// `op` is associative: `op(op(a, b), c) = op(a, op(b, c))` for all `a`, `b`, `c`.
pub trait Semigroup {
    type Value;
    fn op(&self, a: &Self::Value, b: &Self::Value) -> Self::Value;
}

/// A monoid on `Value`.
///
/// # Contract
/// `id()` is the identity: `op(id(), a) = op(a, id()) = a` for all `a`.
pub trait Monoid: Semigroup {
    fn id(&self) -> Self::Value;
}

/// A group on `Value`.
///
/// # Contract
/// `inv(a)` is the inverse of `a`: `op(a, inv(a)) = op(inv(a), a) = id()` for all `a`.
pub trait Group: Monoid {
    fn inv(&self, a: &Self::Value) -> Self::Value;
}

/// A map `act: U x T -> T`, applying the operator `f` in `U` to `x` in `T`.
///
/// # Contract
/// None here. The structures that use an action state the laws they require.
pub trait Action<T, U> {
    fn act(&self, f: &U, x: &T) -> T;
}

/// A semiring on `Value`.
///
/// # Contract
/// - `(Value, zero, add)` is a commutative monoid, and `(Value, one, mul)` is a monoid.
/// - `mul` distributes over `add`: `mul(add(a, b), c) = add(mul(a, c), mul(b, c))` and
///   `mul(a, add(b, c)) = add(mul(a, b), mul(a, c))` for all `a`, `b`, `c`.
/// - `zero()` is abosorbing: `mul(zero(), a) = mul(a, zero()) = zero()` for all `a`.
pub trait Semiring {
    type Value;
    fn zero(&self) -> Self::Value;
    fn one(&self) -> Self::Value;
    fn add(&self, a: &Self::Value, b: &Self::Value) -> Self::Value;
    fn mul(&self, a: &Self::Value, b: &Self::Value) -> Self::Value;
}

/// A ring on `Value`.
///
/// # Contract
/// `neg(a)` is the additive inverse of `a`: add(a, neg(a)) = zero()` for all `a`.
pub trait Ring: Semiring {
    fn neg(&self, a: &Self::Value) -> Self::Value;
}

/// A skew field on `Value`.
///
/// # Contract
/// `zero() != one()`, and `inv(a)` is the multiplicative inverse of every `a != zero()`:
/// `mul(a, inv(a)) = one()`. `inv(zero())` is not defined; implementations may panic.
pub trait SkewField: Ring {
    fn inv(&self, a: &Self::Value) -> Self::Value;
}

/// A field on `Value`, a skew field whose `mul` is commutative.
pub trait Field: SkewField + Commutative {}
impl<R: SkewField + Commutative> Field for R {}

/// A field with a choice of primitive roots of unity.
///
/// # Contract
/// `root_of_unity(n)` is `Some(ω)` for a primitive `n`-th root of unity `ω` if the field has one,
/// and`None` otherwise.
pub trait RootOfUnity: Field {
    fn root_of_unity(&self, n: usize) -> Option<Self::Value>;
}

/// A marker for commutativity.
///
/// # Contract
/// For a semigroup, `op(a, b) = op(b, a)` for all `a`, `b`. For a semiring, whose `add` is always
/// commutative, `mul(a, b) = mul(b, a)` for all `a`, `b`.
pub trait Commutative {}

/// A marker for idempotency.
///
/// # Contract
/// For a semigroup, `op(a, a) = a` for all `a`.
pub trait Idempotent {}

/// A type with a distinguished element `zero`.
///
/// # Contract
/// `zero` is the additive identity of `T`, i.e. `T::zero() + a = a + T::zero() = a` for all `a`.
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

/// A type with a least and a greatest element.
///
/// # Contract
/// `min_value` and `max_value` are the least and the greatest element of `T` under `Ord`,
/// i.e. `T::min_value() <= a <= T::max_value()` for all `a`.
pub trait Bounded: Ord {
    fn min_value() -> Self;
    fn max_value() -> Self;
}

/// A type whose nonzero elements have multiplicative inverses.
///
/// # Contract
/// `inv` returns the multiplicative inverse in `T`, i.e. `a * inv(a) = inv(a) * a = one` for all
/// `a != zero`.
/// `inv(zero())` is not defined; implementations may panic.
pub trait Inv {
    type Output;
    fn inv(&self) -> Self::Output;
}
