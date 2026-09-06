use crate::algebra::lcm::Lcm;
use crate::algebra::{Ring, Semiring};
use crate::convolution::{InverseTransform, Transform};
use crate::num::primes::primes;

impl<R: Semiring> Transform<R> for Lcm<usize> {
    /// The divisor zeta transform `(Zf)(d) = Σ_{x: x|d} f(x)`.
    ///
    /// # Definition
    /// The rows of `Z` are the characters `x -> [x | d]` of `(N, 1, lcm)`:
    /// `lcm(x, y) | d` iff `x | d` and `y | d`, which is `Z(f * g) = Z(f) .* Z(g)`.
    /// They take values in `{0, 1}`, so they exist in every semiring; hence the bound `Semiring`.
    /// `Z` restricted to `[1, n)` is the Kronecker product of the zeta transforms of the chains
    /// `(N, max)` over the primes.
    ///
    /// Index `0` is the top of the divisibility order: `x | 0` for every `x`,
    /// and `0 | d` iff `d = 0`, so `(Zf)(0) = Σ_x f(x)` and `f(0)` contributes to no other
    /// component.
    ///
    /// # Contract
    /// `[0, n)` is not closed under `lcm`. `convolve` returns the product of `R[Lcm]` projected to
    /// `[0, n)`: the terms with `lcm(x, y) >= n` are dropped.
    ///
    /// # Complexity
    /// - Time: O(n log log n)
    /// - Space: O(n)
    fn transform(&self, ring: &R, f: &mut [R::Value]) {
        let n = f.len();
        if n > 1 {
            let (zero, rest) = f.split_first_mut().unwrap();
            for x in rest.iter() {
                *zero = ring.add(zero, x);
            }
        }
        for p in primes(n) {
            for i in 1..n.div_ceil(p) {
                f[i * p] = ring.add(&f[i * p], &f[i]);
            }
        }
    }
}

impl<R: Ring> InverseTransform<R> for Lcm<usize> {
    /// The divisor Möbius transform `Z^{-1}`, `(Z^{-1} g)(d) = Σ_{x: x|d} μ(d / x) g(x)`.
    ///
    /// # Definition
    /// `Z` is unitriangular, so `Z^{-1}` exists over every ring; its entries are the function `μ`.
    /// It needs subtraction but no division, hence the bound `Ring`.
    ///
    /// # Complexity
    /// - Time: O(n log log n)
    /// - Space: O(n)
    fn inverse_transform(&self, ring: &R, f: &mut [R::Value]) {
        let n = f.len();
        for p in primes(n) {
            for i in (1..n.div_ceil(p)).rev() {
                f[i * p] = ring.add(&f[i * p], &ring.neg(&f[i]));
            }
        }
        if n > 1 {
            let (zero, rest) = f.split_first_mut().unwrap();
            for x in rest.iter() {
                *zero = ring.add(zero, &ring.neg(x));
            }
        }
    }
}
