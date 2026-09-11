use crate::algebra::{Commutative, Group, Monoid, Ring, Semigroup, Semiring};

/// The additive group of a semiring `R`, i.e. `R` with its multiplication forgotten.
///
/// # Definition
/// `(R::Value, zero, add)` is a commutative monoid by the definition of a semiring,
/// and `(R::Value, zero, add, neg)` is a commutative group when `R` is a ring.
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
