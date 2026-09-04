use crate::algebra::{Monoid, Semiring};

/// The multiplicative monoid of a semiring `R`, i.e. `R` with its addition forgotten.
///
/// # Definition
/// `(R::Value, one, mul)` is a monoid by the definition of a semiring,
/// and `(R::Value, one, mul)` is a group when `R` is a ring.
#[derive(Clone, Copy, Default)]
pub struct Multiplicative<R>(pub R);

impl<R: Semiring> Monoid for Multiplicative<R> {
    type Value = R::Value;
    fn id(&self) -> R::Value {
        self.0.one()
    }
    fn op(&self, a: &R::Value, b: &R::Value) -> R::Value {
        self.0.mul(a, b)
    }
}
