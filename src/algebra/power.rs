use crate::algebra::{Monoid, Semiring};

/// Raises `x` to the power of `exp` in a monoid by repeated squaring.
///
/// # Definition
/// `x^0 = id` and `x^n = x^{n-1} * x`.
///
/// # Complexity
/// - Time: O(log exp)
/// - Space: O(1)
pub fn pow<M: Monoid>(monoid: &M, x: &M::Value, exp: u64) -> M::Value {
    let mut y = monoid.id();
    for i in (0..u64::BITS - exp.leading_zeros()).rev() {
        y = monoid.op(&y, &y);
        if exp >> i & 1 == 1 {
            y = monoid.op(&y, x);
        }
    }
    y
}

/// The geometric sum `1 + x + ... + x^{n-1}` in a semiring.
///
/// # Definition
/// `Σ_{k=0}^{n-1} x^k`, the empty sum being `zero`.
///
/// # Complexity
/// - Time: O(log n)
/// - Space: O(1)
pub fn geometric_sum<R: Semiring>(ring: &R, x: &R::Value, n: u64) -> R::Value {
    let mut y = ring.zero();
    let mut pow = ring.one();
    for i in (0..u64::BITS - n.leading_zeros()).rev() {
        (y, pow) = (ring.add(&y, &ring.mul(&pow, &y)), ring.mul(&pow, &pow));
        if n >> i & 1 == 1 {
            (y, pow) = (ring.add(&ring.one(), &ring.mul(x, &y)), ring.mul(&pow, x));
        }
    }
    y
}
