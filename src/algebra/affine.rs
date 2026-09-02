use crate::algebra::{Action, Monoid, One, Zero};

/// The affine monoid on `T`, given by its own `std::ops::Add, std::ops::Mul`.
pub struct Affine<T>(std::marker::PhantomData<T>);
impl<T> Affine<T> {
    /// Creates an affine structure.
    pub const fn new() -> Self {
        Self(std::marker::PhantomData)
    }
}
impl<T> Default for Affine<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Clone for Affine<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for Affine<T> {}

impl<T: Clone + std::ops::Add<Output = T> + std::ops::Mul<Output = T> + Zero + One> Monoid
    for Affine<T>
{
    /// `(a, b)` means `x -> ax + b`.
    type Value = (T, T);
    fn id(&self) -> (T, T) {
        (T::one(), T::zero())
    }
    fn op(&self, a: &(T, T), b: &(T, T)) -> (T, T) {
        (
            b.0.clone() * a.0.clone(),
            b.0.clone() * a.1.clone() + b.1.clone(),
        )
    }
}

impl<T: Clone + std::ops::Add<Output = T> + std::ops::Mul<Output = T> + Zero + One>
    Action<T, (T, T)> for Affine<T>
{
    fn act(&self, f: &(T, T), x: &T) -> T {
        f.0.clone() * x.clone() + f.1.clone()
    }
}
