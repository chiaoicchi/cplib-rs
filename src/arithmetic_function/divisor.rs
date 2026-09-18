use crate::algebra::{Ring, Semiring};
use crate::num::prime::primes;

/// The divisor zeta transform of `f`, in place.
///
/// # Definition
/// `(Zf)(d) = Σ_{x|d} f(x)`. `x | 0` for every `x`, so `(Zf)(0) = Σ_x f(x)`, and `0 | d` only for
/// `d = 0`, so `f(0)` contributes to no other component.
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
/// `(Z^{-1}f)(d) = Σ_{x|d} μ(d/x) g(x)`, where `μ` is the Mobius function. `Z` is unitriangular in
/// the divisibility order, so `Z^{-1}` exists over every ring; hence the bound `Ring`.
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

/// the lcm convolution of `f` and `g`.
///
/// # Definition
/// `(fg)(d) = Σ_{lcm(x,y)=d} f(x) g(y)` for `d < n`, the product of monoid algebra `R[(N, lcm)]`
/// with the terms of `lcm(x, y) >= n` discard; `lcm(x, 0) = 0`, so `0` absorbs. Since
/// `lcm(x, y) | d` iff `x | d` and `y | d`, [`divisor_zeta`] carries it to the pointwise product:
/// `Z(fg) = Z(f) Z(g)`.
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
