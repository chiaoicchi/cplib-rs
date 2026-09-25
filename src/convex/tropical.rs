//! Tropical semirings, and the products of convex polynomials over them.
//!
//! A tropical semiring is a semiring `R` whose `add` is selective, `add(a, b)` being `a` or `b` for
//! all `a`, `b`, so that `a <= b` iff `add(a, b) = a` is a total order, and whose `mul` is
//! commutative and cancellative in that order: `a <= b` iff `mul(a, c) <= mul(b, c)` for all `a`,
//! `b` and all `c != zero()`. A polynomial `f` over `R` is stored with `f(i)`, the coefficient of
//! `x^i`, at index `i`, and is convex if `mul(f(i + 1), f(i + 1)) <= mul(f(i), f(i + 2))` for every
//! `i`. Throughout this module `n` and `m` are the lengths of the operands `f` and `g`.

use std::marker::PhantomData;

use crate::algebra::{Bounded, Commutative, Idempotent, Semiring, Zero};
use crate::convex::monge::monotone_minima;

/// The tropical semiring of `T` under `min` and `+`.
///
/// # Definition
/// `∞` is `T::max_value()`. `add(a, b) = min(a, b)`, `zero()` is `∞`, `mul(a, b) = a + b` if
/// neither is `∞` and `∞` otherwise, and `one()` is `T::zero()`.
///
/// # Contract
/// - `a + b` neither overflows nor equals `∞` for all `a`, `b` other than `∞` that are multiplied.
pub struct MinPlus<T>(PhantomData<T>);
impl<T> MinPlus<T> {
    /// The tropical semiring of `T`.
    pub const fn new() -> Self {
        Self(PhantomData)
    }
}
impl<T> Default for MinPlus<T> {
    fn default() -> Self {
        Self::new()
    }
}
impl<T> Clone for MinPlus<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for MinPlus<T> {}

impl<T: Copy + Ord + Bounded + Zero + std::ops::Add<Output = T>> Semiring for MinPlus<T> {
    type Value = T;
    fn zero(&self) -> T {
        T::max_value()
    }
    fn one(&self) -> T {
        T::zero()
    }
    fn add(&self, a: &T, b: &T) -> T {
        *a.min(b)
    }
    fn mul(&self, a: &T, b: &T) -> T {
        if *a == T::max_value() || *b == T::max_value() {
            T::max_value()
        } else {
            *a + *b
        }
    }
}
impl<T> Commutative for MinPlus<T> {}
impl<T> Idempotent for MinPlus<T> {}

/// The tropical semiring of `T` under `max` and `+`.
///
/// # Definition
/// `-∞` is `T::min_value()`. `add(a, b) = max(a, b)`, `zero()` is `-∞`, `mul(a, b) = a + b` if
/// neither is `-∞` and `-∞` otherwise, and `one()` is `T::zero()`.
///
/// # Contract
/// - `a + b` neither overflows nor equals `-∞` for all `a`, `b` other than `-∞` that are
/// multiplied.
pub struct MaxPlus<T>(PhantomData<T>);
impl<T> MaxPlus<T> {
    /// The tropical semiring of `T`.
    pub const fn new() -> Self {
        Self(PhantomData)
    }
}
impl<T> Default for MaxPlus<T> {
    fn default() -> Self {
        Self::new()
    }
}
impl<T> Clone for MaxPlus<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for MaxPlus<T> {}

impl<T: Copy + Ord + Bounded + Zero + std::ops::Add<Output = T>> Semiring for MaxPlus<T> {
    type Value = T;
    fn zero(&self) -> Self::Value {
        T::min_value()
    }
    fn one(&self) -> Self::Value {
        T::zero()
    }
    fn add(&self, a: &T, b: &T) -> T {
        *a.max(b)
    }
    fn mul(&self, a: &T, b: &T) -> T {
        if *a == T::min_value() || *b == T::min_value() {
            T::min_value()
        } else {
            *a + *b
        }
    }
}
impl<T> Commutative for MaxPlus<T> {}
impl<T> Idempotent for MaxPlus<T> {}

/// The product of convex `f` and `g` in `R[x]`.
///
/// # Definition
/// `(fg)(k) = Σ_{i+j=k} f(i) g(j)`, with `add` as the sum and `mul` as the product. The result has
/// length `n + m - 1`, and is empty if `n = 0` or `m = 0`.
///
/// # Contract
/// - `R` is a tropical semiring.
///
/// # Complexity
/// - Time: O(n + m)
/// - Space: O(n + m)
///
/// # Panics
/// Panics if `f` or `g` is not convex.
/// Panics if some element of `f` or `g` is `zero()`.
pub fn tropical_convolve_convex_convex<
    R: Semiring<Value: Clone + PartialEq> + Commutative + Idempotent,
>(
    ring: &R,
    f: &[R::Value],
    g: &[R::Value],
) -> Vec<R::Value> {
    let le = |a: &R::Value, b: &R::Value| ring.add(a, b) == *a;
    let zero = ring.zero();
    assert!(
        f.iter().chain(g).all(|x| *x != zero),
        "f and g must not contain zero()"
    );
    assert!(
        f.windows(3)
            .all(|w| le(&ring.mul(&w[1], &w[1]), &ring.mul(&w[0], &w[2]))),
        "f must be convex"
    );
    assert!(
        g.windows(3)
            .all(|w| le(&ring.mul(&w[1], &w[1]), &ring.mul(&w[0], &w[2]))),
        "g must be convex"
    );
    let (n, m) = (f.len(), g.len());
    if n == 0 || m == 0 {
        return vec![];
    }
    let (mut i, mut j) = (0, 0);
    let mut h = Vec::with_capacity(n + m - 1);
    h.push(ring.mul(&f[0], &g[0]));
    while i + 1 < n || j + 1 < m {
        if j + 1 == m || (i + 1 < n && le(&ring.mul(&f[i + 1], &g[j]), &ring.mul(&f[i], &g[j + 1])))
        {
            i += 1;
        } else {
            j += 1;
        }
        h.push(ring.mul(&f[i], &g[j]));
    }
    h
}

/// The product of convex `f` and arbitrary `g` in `R[x]`.
///
/// # Definition
/// `(fg)(k) = Σ_{i+j=k} f(i) g(j)`, with `add` as the sum and `mul` as the product. The result has
/// length `n + m - 1`, and is empty if `n = 0` or `m = 0`.
///
/// # Contract
/// - `R` is a tropical semiring.
///
/// # Complexity
/// - Time: O((n + m) log (n + m))
/// - Space: O(n + m)
///
/// # Panics
/// Panics if `f` is not convex.
/// Panics if some element of `f` or `g` is `zero()`.
pub fn tropical_convolve_convex_arbitrary<
    R: Semiring<Value: Clone + PartialEq> + Commutative + Idempotent,
>(
    ring: &R,
    f: &[R::Value],
    g: &[R::Value],
) -> Vec<R::Value> {
    struct Entry<'a, R: Semiring> {
        ring: &'a R,
        value: R::Value,
    }
    impl<R: Semiring<Value: PartialEq>> PartialEq for Entry<'_, R> {
        fn eq(&self, other: &Self) -> bool {
            self.value == other.value
        }
    }
    impl<R: Semiring<Value: PartialEq>> PartialOrd for Entry<'_, R> {
        fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
            use std::cmp::Ordering::{Equal, Greater, Less};
            Some(if self.value == other.value {
                Equal
            } else if self.ring.add(&self.value, &other.value) == self.value {
                Less
            } else {
                Greater
            })
        }
    }
    let le = |a: &R::Value, b: &R::Value| ring.add(a, b) == *a;
    let zero = ring.zero();
    assert!(
        f.iter().chain(g).all(|x| *x != zero),
        "f and g must not contain zero()"
    );
    assert!(
        f.windows(3)
            .all(|w| le(&ring.mul(&w[1], &w[1]), &ring.mul(&w[0], &w[2]))),
        "f must be convex"
    );
    let (n, m) = (f.len(), g.len());
    if n == 0 || m == 0 {
        return vec![];
    }
    let argmin = monotone_minima(n + m - 1, m, |k, j| Entry {
        ring,
        value: if j <= k && k - j < n {
            ring.mul(&f[k - j], &g[j])
        } else {
            ring.zero()
        },
    });
    argmin
        .iter()
        .enumerate()
        .map(|(k, &j)| ring.mul(&f[k - j], &g[j]))
        .collect()
}
