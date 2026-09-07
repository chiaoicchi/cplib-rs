use crate::algebra::Semiring;
use crate::algebra::canonical::Canonical;
use crate::algebra::multiplicative::Multiplicative;
use crate::convolution::Convolution;

impl<R: Semiring> Convolution<R> for Multiplicative<Canonical<usize>> {
    /// The Dirichlet convolution `(f * g)(m) = Σ_{ab=m} f(a)g(b)`, tabulated on `[0, n)`.
    ///
    /// No transform diagonalizes `R[(N, x)]` truncated to `[0, n)`: it is
    /// `R[x_p: p prime]` modulo the monomials of index `>= n`, in which every `x_p` is
    /// nilpotent, so the algebra is not semisimple. The product is computed directly.
    ///
    /// # Complexity
    /// - Time: O(n log n)
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `f.len() != g.len()`.
    fn convolve(&self, ring: &R, f: &[R::Value], g: &[R::Value]) -> Vec<R::Value> {
        assert_eq!(
            f.len(),
            g.len(),
            "length must agree: {} != {}",
            f.len(),
            g.len()
        );
        let n = f.len();
        let mut h: Vec<R::Value> = (0..n).map(|_| ring.zero()).collect();
        for a in 1..n {
            for b in 1..n.div_ceil(a) {
                h[a * b] = ring.add(&h[a * b], &ring.mul(&f[a], &g[b]));
            }
        }
        h
    }
}
