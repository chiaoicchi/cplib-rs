use crate::algebra::{Ring, Semiring};
use crate::num::prime::primes;

/// The multiple zeta transform of `f`, in place.
///
/// # Definition
/// `(Zf)(d) = Σ_{d|x} f(x)`. `d | 0` for every `d`, so `f(0)` contributes to every component, and
/// `0 | x` only for `x = 0` so `(Zf)(0) = f(0)`.
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
/// `(Z^{-1}g)(d) = Σ_{d|x} μ(x/d) g(x)`, where `μ` is the Mobius function. `Z` is unitriangular in
/// the divisibility order, so `Z^{-1}` exists over every ring; hence the bound `Ring`.
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

/// the gcd convolution of `f` and `g`.
///
/// # Definition
/// `(fg)(d) = Σ_{gcd(x,y)=d} f(x) g(y)`, the product of the monoid algebra `R[(N, gcd)]`, whose
/// identity is `0`. Since `d | x` and `d | y` iff `d | gcd(x, y)`, [`multiple_zeta`] carries it to
/// the pointwise product: `Z(fg) = Z(f) Z(g)`.
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
