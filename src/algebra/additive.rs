use crate::algebra::{Commutative, Group, Monoid};

/// A additive group of `u64`, i.e. the integers modulo 2^64.
///
/// The operation is wrapping addition, the identity is `0`, and the inverse is wrapping negation.
/// The operation is commutative.
pub struct Additive;
impl Monoid for Additive {
    type Value = u64;
    fn id(&self) -> u64 {
        0
    }
    fn op(&self, a: &u64, b: &u64) -> u64 {
        a.wrapping_add(*b)
    }
}
impl Group for Additive {
    fn inv(&self, a: &u64) -> u64 {
        a.wrapping_neg()
    }
}
impl Commutative for Additive {}
