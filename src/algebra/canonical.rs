//! Structures given by the operators of a type, and the monoids of a semiring.

use crate::algebra::{
    Action, Commutative, Group, Inv, Monoid, One, Ring, Semigroup, Semiring, SkewField, Zero,
};

/// The structures of `T` given by its own operators.
///
/// # Definition
/// `zero()`, `one()`, `add`, `mul`, `neg` and `inv` are `T::zero()`, `T::one()`, `+`, `*`, unary
/// `-` and [`Inv::inv`], for each of them that `T` has.
///
/// # Contract
/// The operators of `T` satisfy the laws of the traits `Canonical<T>` implements.
pub struct Canonical<T>(std::marker::PhantomData<T>);
impl<T> Canonical<T> {
    /// The structure of `T` given by its own operators.
    pub const fn new() -> Self {
        Self(std::marker::PhantomData)
    }
}
impl<T> Default for Canonical<T> {
    fn default() -> Self {
        Self::new()
    }
}
impl<T> Clone for Canonical<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for Canonical<T> {}

impl<T: Clone + std::ops::Add<Output = T> + std::ops::Mul<Output = T> + Zero + One> Semiring
    for Canonical<T>
{
    type Value = T;
    fn zero(&self) -> T {
        T::zero()
    }
    fn one(&self) -> T {
        T::one()
    }
    fn add(&self, a: &T, b: &T) -> T {
        a.clone() + b.clone()
    }
    fn mul(&self, a: &T, b: &T) -> T {
        a.clone() * b.clone()
    }
}
impl<
    T: Clone
        + std::ops::Add<Output = T>
        + std::ops::Neg<Output = T>
        + std::ops::Mul<Output = T>
        + Zero
        + One,
> Ring for Canonical<T>
{
    fn neg(&self, a: &T) -> T {
        -a.clone()
    }
}
impl<
    T: Clone
        + std::ops::Add<Output = T>
        + std::ops::Mul<Output = T>
        + std::ops::Neg<Output = T>
        + Zero
        + One
        + Inv<Output = T>,
> SkewField for Canonical<T>
{
    fn inv(&self, a: &T) -> T {
        Inv::inv(a)
    }
}

/// The monoid of a semiring `R` under its addition.
///
/// # Definition
/// `op = add` and `id() = zero()`, with `inv = neg` if `R` is a ring.
#[derive(Clone, Copy, Default)]
pub struct Additive<R>(pub R);

impl<R: Semiring> Semigroup for Additive<R> {
    type Value = R::Value;
    fn op(&self, a: &R::Value, b: &R::Value) -> R::Value {
        self.0.add(a, b)
    }
}
impl<R: Semiring> Monoid for Additive<R> {
    fn id(&self) -> R::Value {
        self.0.zero()
    }
}
impl<R: Ring> Group for Additive<R> {
    fn inv(&self, a: &R::Value) -> R::Value {
        self.0.neg(a)
    }
}
impl<R: Semiring> Commutative for Additive<R> {}

/// The monoid of a semiring `R` under its multiplication.
#[derive(Clone, Copy, Default)]
pub struct Multiplicative<R>(pub R);
impl<R: Semiring> Semigroup for Multiplicative<R> {
    type Value = R::Value;
    fn op(&self, a: &R::Value, b: &R::Value) -> R::Value {
        self.0.mul(a, b)
    }
}
impl<R: Semiring> Monoid for Multiplicative<R> {
    fn id(&self) -> R::Value {
        self.0.one()
    }
}
impl<R: Semiring + Commutative> Commutative for Multiplicative<R> {}

/// The monoid of affine maps on a semiring `R`, under composition.
///
/// # Definition
/// An element `(a, b)` is the map `x -> ax + b` on `R`. `op(f, g) = g .* f`, the composite that
/// applies `f` first, and `id()` is `(1, 0)`. As an action, `act(f, x) = f(x)`.
#[derive(Clone, Copy, Default)]
pub struct Affine<R>(pub R);
impl<R: Semiring> Semigroup for Affine<R> {
    type Value = (R::Value, R::Value);
    fn op(&self, f: &(R::Value, R::Value), g: &(R::Value, R::Value)) -> (R::Value, R::Value) {
        (
            self.0.mul(&g.0, &f.0),
            self.0.add(&self.0.mul(&g.0, &f.1), &g.1),
        )
    }
}
impl<R: Semiring> Monoid for Affine<R> {
    fn id(&self) -> (R::Value, R::Value) {
        (self.0.one(), self.0.zero())
    }
}
impl<R: Semiring> Action<R::Value, (R::Value, R::Value)> for Affine<R> {
    fn act(&self, f: &(R::Value, R::Value), x: &R::Value) -> R::Value {
        self.0.add(&self.0.mul(&f.0, x), &f.1)
    }
}
