use crate::algebra::{Commutative, Group, Monoid, Zero};

/// The additive group of `T`, given by its own `std::ops::Add`.
pub struct Additive<T>(std::marker::PhantomData<T>);
impl<T> Additive<T> {
    /// Creates an additive structure.
    pub const fn new() -> Self {
        Self(std::marker::PhantomData)
    }
}
impl<T> Default for Additive<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Clone for Additive<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for Additive<T> {}

impl<T: Clone + std::ops::Add<Output = T> + Zero> Monoid for Additive<T> {
    type Value = T;
    fn id(&self) -> T {
        T::zero()
    }
    fn op(&self, a: &T, b: &T) -> T {
        a.clone() + b.clone()
    }
}

impl<T: Clone + std::ops::Add<Output = T> + std::ops::Neg<Output = T> + Zero> Group
    for Additive<T>
{
    fn inv(&self, a: &T) -> T {
        -a.clone()
    }
}
impl<T: Clone + std::ops::Add<Output = T>> Commutative for Additive<T> {}
