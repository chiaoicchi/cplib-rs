use crate::algebra::xor::Xor;
use crate::algebra::{Field, Ring};
use crate::convolution::{InverseTransform, Transform};

impl<R: Ring> Transform<R> for Xor<usize> {
    /// The Walsh-Hadamard transform `(Hf)(s) = Σ_x (-1)^{<s, x>} f(x)`.
    ///
    /// # Definition
    /// The rows of `H` are the characters `x -> (-1)^{<s, x>}` of `(F_2)^n`, which exist iff
    /// `R` has an element `-1` with `(-1)^2 = 1` and `-1 != 1`, i.e. iff `R` is a ring of
    /// characteristic other than `2`; hence the bound `Ring`. `H` is the `n`-fold Kronecker product
    /// of `H_1 = [[1, 1], [1, -1]]`.
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
            for i in (0..f.len()).filter(|i| i & w == 0) {
                (f[i], f[i | w]) = (
                    ring.add(&f[i], &f[i | w]),
                    ring.add(&f[i], &ring.neg(&f[i | w])),
                );
            }
            w <<= 1;
        }
    }
}
impl<R: Field> InverseTransform<R> for Xor<usize> {
    /// The inverse Walsh-Hadamard transform `H^{-1} = N^{-1}H`, where `N = 2^n`.
    ///
    /// # Definition
    /// `H^2 = N id`, since `Σ_s (-1)^{<s, x + y>} = N [x = y]`; hence `H^{-1} = N^{-1} H`.
    /// This requires `N = 2^n`, and therefore `2`, to be invertible in `R`, which is the reason for
    /// the bound `Field`.
    ///
    /// # Complexity
    /// - Time: O(2^n n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `f.len()` is not a power of two.
    fn inverse_transform(&self, ring: &R, f: &mut [R::Value]) {
        self.transform(ring, f);
        let half = ring.inv(&ring.add(&ring.one(), &ring.one()));
        let mut c = ring.one();
        for _ in 0..f.len().trailing_zeros() {
            c = ring.mul(&c, &half);
        }
        for x in f.iter_mut() {
            *x = ring.mul(&c, x);
        }
    }
}
