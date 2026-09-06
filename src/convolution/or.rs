use crate::algebra::or::Or;
use crate::algebra::{Ring, Semiring};
use crate::convolution::{InverseTransform, Transform};

impl<R: Semiring> Transform<R> for Or<usize> {
    /// The subset zeta transform `(Zf)(s) = Σ_{x: x|s=s} f(x)`.
    ///
    /// # Definition
    /// The rows of `Z` are the characters `x -> [x | s = s]` of `({0, 1}, 0, |)^n`:
    /// `(x | y) | s = s` iff `x | s = s` and `y | s = s`, which is `Z(f * g) = Z(f) .* Z(g)`.
    /// They take values in `{0, 1}`, so they exist in every semiring; hence the bound `Semiring`.
    /// `Z` is the `n`-fold Kronecker product of `Z_1 = [[1, 0], [1, 1]]`.
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
                    f[j | w] = ring.add(&f[j | w], &f[j]);
                }
            }
            w <<= 1;
        }
    }
}

impl<R: Ring> InverseTransform<R> for Or<usize> {
    /// The subset Möbius transform `Z^{-1}`.
    ///
    /// # Definition
    /// `Z_1^{-1} = [[1, 0], [-1, 1]]`, so `Z^{-1}` is its `n`-fold Kronecker product. It needs
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
                    f[j | w] = ring.add(&f[j | w], &ring.neg(&f[j]));
                }
            }
            w <<= 1;
        }
    }
}
