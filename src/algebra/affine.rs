use crate::algebra::{Action, Monoid, Semiring};

/// The monoid of affine maps over a semiring `R`, under composition.
///
/// # Definition
/// An element `(a, b)` represents the map `x -> a * x + b`.
/// `op(f, g)(x)` is the composite `g(f(x))`, i.e. `f` is applied first.
#[derive(Clone, Copy, Default)]
pub struct Affine<R>(pub R);

impl<R: Semiring> Monoid for Affine<R> {
    type Value = (R::Value, R::Value);
    fn id(&self) -> (R::Value, R::Value) {
        (self.0.one(), self.0.zero())
    }
    fn op(&self, f: &(R::Value, R::Value), g: &(R::Value, R::Value)) -> (R::Value, R::Value) {
        (
            self.0.mul(&g.0, &f.0),
            self.0.add(&self.0.mul(&g.0, &f.1), &g.1),
        )
    }
}

impl<R: Semiring> Action<R::Value, (R::Value, R::Value)> for Affine<R> {
    fn act(&self, f: &(R::Value, R::Value), x: &R::Value) -> R::Value {
        self.0.add(&self.0.mul(&f.0, x), &f.1)
    }
}
