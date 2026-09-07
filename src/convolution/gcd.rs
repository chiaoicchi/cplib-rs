use crate::algebra::gcd::Gcd;
use crate::algebra::{Ring, Semiring};
use crate::convolution::{InverseTransform, Transform};
use crate::num::prime::primes;

impl<R: Semiring> Transform<R> for Gcd<usize> {
    /// The multiple zeta transform `(Zf)(d) = Σ_{x: d|x} f(x)`.
    ///
    /// # Definition
    /// The rows of `Z` are the characters `x -> [d | x]` of `(N, 0, gcd)`:
    /// `d | gcd(x, y)` iff `d | x` and `d | y`, which is `Z(f * g) = Z(f) .* Z(g)`.
    /// They take values in `{0, 1}`, so they exist in every semiring; hence the bound `Semiring`.
    ///
    /// Under `x -> (v_p(x))_p` the positive integers are the direct product of the chains
    /// `(N, min)` over the primes `p`, so `Z` restricted to `[1, n)` is the Kronecker product of
    /// the zeta transforms of these chains.
    ///
    /// Index `0` is the identity of `gcd`: `0 | x` iff `x = 0`, and `d | 0` for every `d`, so
    /// `(Zf)(0) = f(0)` and `(Zf)(d) = f(0) + Σ_{x >= 1: d|x} f(x)` for `d >= 1`.
    ///
    /// # Complexity
    /// - Time: O(n log log n)
    /// - Space: O(n)
    fn transform(&self, ring: &R, f: &mut [R::Value]) {
        let n = f.len();
        for p in primes(n - 1) {
            for i in (1..n.div_ceil(p)).rev() {
                f[i] = ring.add(&f[i], &f[i * p]);
            }
        }
        if n > 1 {
            let (zero, rest) = f.split_first_mut().unwrap();
            for x in rest {
                *x = ring.add(x, zero);
            }
        }
    }
}

impl<R: Ring> InverseTransform<R> for Gcd<usize> {
    /// The multiple Möbius transform `Z^{-1}`, `(Z^{-1} g)(d) = Σ_{x: d|x} μ(x / d) g(x)`.
    ///
    /// # Definition
    /// `Z` is unitriangular, so `Z^{-1}` exists over every ring; its entries are the Möbius
    /// function `μ`, which takes values in `{0, 1, -1}`. It needs subtraction but no division,
    /// hence the bound `Ring` rather than `Field`. As `Z` is the Kronecker product of chain zeta
    /// transforms, `Z^{-1}` is the Kronecker product of their inverses.
    ///
    /// # Complexity
    /// - Time: O(n log log n)
    /// - Space: O(n)
    fn inverse_transform(&self, ring: &R, f: &mut [R::Value]) {
        let n = f.len();
        if n > 1 {
            let (zero, rest) = f.split_first_mut().unwrap();
            let neg = ring.neg(zero);
            for x in rest {
                *x = ring.add(x, &neg);
            }
        }
        for p in primes(n - 1) {
            for i in 1..n.div_ceil(p) {
                f[i] = ring.add(&f[i], &ring.neg(&f[i * p]));
            }
        }
    }
}
