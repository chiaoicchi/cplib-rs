use crate::algebra::{Field, Inv, One, Ring, Semiring, Zero};

/// The canonical structure of `T`, given by its own operators.
pub struct Canonical<T>(std::marker::PhantomData<T>);
impl<T> Canonical<T> {
    /// Creates a Canonical structure.
    pub const fn new() -> Self {
        Self(std::marker::PhantomData)
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
> Field for Canonical<T>
{
    fn inv(&self, a: &T) -> T {
        Inv::inv(a)
    }
}
