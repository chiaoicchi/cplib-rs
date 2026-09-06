use crate::algebra::and::And;
use crate::algebra::{Ring, Semiring};
use crate::convolution::{InverseTransform, Transform};

impl<R: Semiring> Transform<R> for And<usize> {
    /// The superset zeta transform `(Zf)(s) = Σ_{x: s&x=s} f(x)`.
    ///
    /// # Definition
    /// The rows of `Z` are the characters `x -> [s & x = s]` of `({0, 1}, &)^n`:
    /// `s & (x & y) = s` iff `s & x = s` and `s & y = s`, which is `Z(f * g) = Z(f) .* Z(g).
    /// They take values in `{0, 1}`, so they exist in every semiring; hence the bound `Semiring`.
    /// `Z` is the `n`-fold Kronecker product of `Z_1 = [[1, 1], [0, 1]]`, i.e.
    /// `Z_1` applied along each of the `n` axes of `f` viewed as a `2 x ... x 2` array.
    ///
    /// # Complexity
    /// - Time: O(2^n n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `f.len()` is not a power of two.
    fn transform(&self, ring: &R, f: &mut [R::Value]) {
        assert!(
            f.len().is_power_of_two(),
            "length must be a power of two: len={}",
            f.len()
        );
        let mut w = 1;
        while w < f.len() {
            for i in (0..f.len()).step_by(w << 1) {
                for j in i..i + w {
                    f[j] = ring.add(&f[j], &f[j | w]);
                }
            }
            w <<= 1;
        }
    }
}
impl<R: Ring> InverseTransform<R> for And<usize> {
    /// The superset Möbius transform `Z^{-1}`.
    ///
    /// # Definition
    /// `Z_1^{-1} = [[1, -1, [0, 1]]`, so `Z^{-1}` is its `n`-fold Kronecker product. It needs
    /// subtraction but no division, hence the bound `Ring` rather than `Field`.
    ///
    /// # Complexity
    /// - Time: O(2^n n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `f.len()` is not a power of two.
    fn inverse_transform(&self, ring: &R, f: &mut [R::Value]) {
        assert!(
            f.len().is_power_of_two(),
            "length must be a power of two: len={}",
            f.len()
        );
        let mut w = 1;
        while w < f.len() {
            for i in (0..f.len()).step_by(w << 1) {
                for j in i..i + w {
                    f[j] = ring.add(&f[j], &ring.neg(&f[j | w]));
                }
            }
            w <<= 1;
        }
    }
}
