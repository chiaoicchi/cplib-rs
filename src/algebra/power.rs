use crate::algebra::{Monoid, Semiring};

/// The power `x^exp` in a monoid.
///
/// # Definition
/// `x^0 = id()` and `x^e = op(x^{e-1}, x)` for `e >= 1`.
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
/// `Σ_{k=0}^{n-1} x^k`, which is `zero()` for `n = 0`.
///
/// # Complexity
/// - Time: O(log n)
/// - Space: O(1)
pub fn geometric_sum<R: Semiring>(ring: &R, x: &R::Value, n: u64) -> R::Value {
    let one = ring.one();
    let mut y = ring.zero();
    let mut pow = ring.one();
    for i in (0..u64::BITS - n.leading_zeros()).rev() {
        (y, pow) = (ring.add(&y, &ring.mul(&pow, &y)), ring.mul(&pow, &pow));
        if n >> i & 1 == 1 {
            (y, pow) = (ring.add(&one, &ring.mul(x, &y)), ring.mul(&pow, x));
        }
    }
    y
}
