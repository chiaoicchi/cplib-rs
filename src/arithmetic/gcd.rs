use crate::algebra::{Commutative, Idempotent, Monoid, Ring, Semigroup, Semiring, Zero};
use crate::arithmetic::prime::primes;

/// The greatest common divisor of `a` and `b`.
///
/// # Definition
/// The common divisor of `a` and `b` that every common divisor divides, with `gcd(a, 0) = a`. In
/// particular `gcd(0, 0) = 0`.
///
/// # Contract
/// `%` is the Euclidean remainder on `T`, and `a`, `b` are non-negative.
///
/// # Complexity
/// - Time: O(log min(a, b))
/// - Space: O(1)
pub fn gcd<T: Clone + PartialEq + Zero + std::ops::Rem<Output = T>>(a: T, b: T) -> T {
    let (mut a, mut b) = (a, b);
    while b != T::zero() {
        (a, b) = (b.clone(), a % b);
    }
    a
}

/// The monoid of `T` under [`gcd`].
///
/// # Definition
/// `(N, 0, gcd)`, the meet of the divisibility order on `N`: `op(a, b) = gcd(a, b)`, and `id()` is
/// `0`.
///
/// # Contract
/// `%` is the Euclidean remainder on `T`, and the values are non-negative.
pub struct Gcd<T>(std::marker::PhantomData<T>);
impl<T> Gcd<T> {
    pub const fn new() -> Self {
        Self(std::marker::PhantomData)
    }
}
impl<T> Default for Gcd<T> {
    fn default() -> Self {
        Self::new()
    }
}
impl<T> Clone for Gcd<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for Gcd<T> {}

impl<T: Clone + PartialEq + std::ops::Rem<Output = T> + Zero> Semigroup for Gcd<T> {
    type Value = T;
    fn op(&self, a: &T, b: &T) -> T {
        gcd(a.clone(), b.clone())
    }
}
impl<T: Clone + PartialEq + std::ops::Rem<Output = T> + Zero> Monoid for Gcd<T> {
    fn id(&self) -> T {
        T::zero()
    }
}
impl<T: Clone + PartialEq + std::ops::Rem<Output = T> + Zero> Commutative for Gcd<T> {}
impl<T: Clone + PartialEq + std::ops::Rem<Output = T> + Zero> Idempotent for Gcd<T> {}

/// The multiple zeta transform of `f`, in place.
///
/// # Definition
/// `(Zf)(d) = Σ_{d|x} f(x)`, where every `d` divides `0` and `0` divides only `0`.
///
/// # Complexity
/// - Time: O(n log log n)
/// - Space: O(n)
pub fn multiple_zeta<R: Semiring>(ring: &R, f: &mut [R::Value]) {
    let n = f.len();
    for p in primes(n.saturating_sub(1)) {
        for i in (1..=(n - 1) / p).rev() {
            f[i] = ring.add(&f[i], &f[i * p]);
        }
    }
    if let Some((zero, rest)) = f.split_first_mut() {
        for x in rest {
            *x = ring.add(x, zero);
        }
    }
}

/// The multiple Mobius transform of `f`, in place, inverting [`multiple_zeta`].
///
/// # Definition
/// `(Z^{-1}f)(0) = f(0)`, and `(Z^{-1}f)(d) = Σ_{d|x, x>=1} μ(x/d) (f(x) - f(0))` for `d >= 1`,
/// where `μ` is the Mobius function.
///
/// # Complexity
/// - Time: O(n log log n)
/// - Space: O(n)
pub fn multiple_mobius<R: Ring>(ring: &R, f: &mut [R::Value]) {
    let Some((zero, rest)) = f.split_first_mut() else {
        return;
    };
    let neg = ring.neg(zero);
    for x in rest {
        *x = ring.add(x, &neg);
    }
    let n = f.len();
    for p in primes(n - 1) {
        for i in 1..=(n - 1) / p {
            f[i] = ring.add(&f[i], &ring.neg(&f[i * p]));
        }
    }
}

/// The gcd convolution of `f` and `g`.
///
/// # Definition
/// `(fg)(d) = Σ_{gcd(x,y)=d} f(x) g(y)`, the product of the monoid algebra of `[0, n)` under `gcd`.
///
/// # Complexity
/// - Time: O(n log log n)
/// - Space: O(n)
///
/// # Panics
/// Panics if `f` and `g` differ in length.
pub fn gcd_convolve<R: Ring>(
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
    multiple_zeta(ring, &mut f);
    multiple_zeta(ring, &mut g);
    for (fi, gi) in f.iter_mut().zip(g.iter()) {
        *fi = ring.mul(fi, gi);
    }
    multiple_mobius(ring, &mut f);
    f
}
