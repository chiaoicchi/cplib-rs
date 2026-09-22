use crate::algebra::{Commutative, Idempotent, Monoid, One, Ring, Semigroup, Semiring, Zero};
use crate::arithmetic::gcd::gcd;
use crate::arithmetic::prime::primes;

/// The least common multiple of `a` and `b`.
///
/// # Definition
/// The common multiple of `a` and `b` that divides every common multiple, with `lcm(a, 0) = 0`.
///
/// # Contract
/// `%` is the Euclidean remainder on `T`, `a`, `b` are non-negative, and `lcm(a, b)` fits in `T`.
///
/// # Complexity
/// - Time: O(log min(a, b))
/// - Space: O(1)
pub fn lcm<
    T: Clone
        + PartialEq
        + Zero
        + std::ops::Mul<Output = T>
        + std::ops::Div<Output = T>
        + std::ops::Rem<Output = T>,
>(
    a: T,
    b: T,
) -> T {
    if a == T::zero() {
        return a;
    }
    a.clone() / gcd(a, b.clone()) * b
}

/// The monoid of `T` under [`lcm`].
///
/// # Definition
/// `(N, 1, lcm)`, the join of the divisibility order on `N`: `op(a, b) = lcm(a, b)`, and `id()` is
/// `1`.
///
/// # Contract
/// `%` is the Euclidean remainder on `T`, the values are non-negative, and every `lcm` fits in `T`.
pub struct Lcm<T>(std::marker::PhantomData<T>);
impl<T> Lcm<T> {
    /// The monoid of `T` under [`lcm`].
    pub const fn new() -> Self {
        Self(std::marker::PhantomData)
    }
}
impl<T> Default for Lcm<T> {
    fn default() -> Self {
        Self::new()
    }
}
impl<T> Clone for Lcm<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for Lcm<T> {}

impl<
    T: Clone
        + PartialEq
        + std::ops::Mul<Output = T>
        + std::ops::Div<Output = T>
        + std::ops::Rem<Output = T>
        + Zero
        + One,
> Semigroup for Lcm<T>
{
    type Value = T;
    fn op(&self, a: &T, b: &T) -> T {
        lcm(a.clone(), b.clone())
    }
}
impl<
    T: Clone
        + PartialEq
        + std::ops::Mul<Output = T>
        + std::ops::Div<Output = T>
        + std::ops::Rem<Output = T>
        + Zero
        + One,
> Monoid for Lcm<T>
{
    fn id(&self) -> T {
        T::one()
    }
}
impl<
    T: Clone
        + PartialEq
        + std::ops::Mul<Output = T>
        + std::ops::Div<Output = T>
        + std::ops::Rem<Output = T>
        + Zero
        + One,
> Commutative for Lcm<T>
{
}
impl<
    T: Clone
        + PartialEq
        + std::ops::Mul<Output = T>
        + std::ops::Div<Output = T>
        + std::ops::Rem<Output = T>
        + Zero
        + One,
> Idempotent for Lcm<T>
{
}

/// The divisor zeta transform of `f`, in place.
///
/// # Definition
/// `(Zf)(d) = Σ_{x|d} f(x)`, where every `x` divides `0` and `0` divides only `0`.
///
/// # Complexity
/// - Time: O(n log log n)
/// - Space: O(n)
pub fn divisor_zeta<R: Semiring>(ring: &R, f: &mut [R::Value]) {
    let Some((zero, rest)) = f.split_first_mut() else {
        return;
    };
    for x in rest.iter() {
        *zero = ring.add(zero, x);
    }
    let n = f.len();
    for p in primes(n - 1) {
        for i in 1..=(n - 1) / p {
            f[i * p] = ring.add(&f[i * p], &f[i]);
        }
    }
}

/// The divisor Mobius transform of `f`, in place, inverting [`divisor_zeta`].
///
/// # Definition
/// `(Z^{-1}f)(d) = Σ_{x|d} μ(d/x) f(x)` for `d >= 1`, where `μ` is the Mobius function, and
/// `(Z^{-1}f)(0) = f(0) - Σ_{d=1,...,n-1} (Z^{-1}f)(d)`.
///
/// # Complexity
/// - Time: O(n log log n)
/// - Space: O(n)
pub fn divisor_mobius<R: Ring>(ring: &R, f: &mut [R::Value]) {
    let n = f.len();
    for p in primes(n.saturating_sub(1)) {
        for i in (1..=(n - 1) / p).rev() {
            f[i * p] = ring.add(&f[i * p], &ring.neg(&f[i]));
        }
    }
    if let Some((zero, rest)) = f.split_first_mut() {
        for x in rest.iter() {
            *zero = ring.add(zero, &ring.neg(x));
        }
    }
}

/// The lcm convolution of `f` and `g`.
///
/// # Definition
/// `(fg)(d) = Σ_{lcm(x,y)=d} f(x) g(y)` for `1 <= d < n`, and `(fg)(0)` collects the terms with
/// `lcm(x, y) = 0` or `lcm(x, y) >= n`: the product of the monoid algebra `R[(N, lcm)/I]`, where
/// `I = {0} ∪ [n, ∞)` is an ideal of `(N, lcm)` and `0` its class.
///
/// # Complexity
/// - Time: O(n log log n)
/// - Space: O(n)
///
/// # Panics
/// Panics if `f` and `g` differ in length.
pub fn lcm_convolve<R: Ring>(
    ring: &R,
    mut f: Vec<R::Value>,
    mut g: Vec<R::Value>,
) -> Vec<R::Value> {
    assert_eq!(
        f.len(),
        g.len(),
        "f and g must have the same length: f={}, g={}",
        f.len(),
        g.len()
    );
    divisor_zeta(ring, &mut f);
    divisor_zeta(ring, &mut g);
    for (fi, gi) in f.iter_mut().zip(g.iter()) {
        *fi = ring.mul(fi, gi);
    }
    divisor_mobius(ring, &mut f);
    f
}
