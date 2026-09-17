use crate::algebra::RootOfUnity;
use crate::algebra::additive::Additive;
use crate::algebra::canonical::Canonical;
use crate::algebra::cyclic::Cyclic;
use crate::convolution::{Convolution, OnlineConvolution, convolve_by_transform};
use crate::divide_and_conquer::cdq;

impl<R: RootOfUnity> Convolution<R> for Additive<Canonical<usize>> {
    /// The polynomial product `(f * g)(k) = Σ_{i + j = k} f(i) g(j)`, tablated on `[0, n)`.
    ///
    /// `R[(N, +)] = R[t]` is computed by the cyclic convolution of a length at least `2n - 1`,
    /// in which the product has no wrap-around (see [`convolve_poly`]).
    ///
    /// # Complexity
    /// - Time: O(n log n)
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `R` has no primitive root of unity of order at least `2n`.
    fn convolve(&self, ring: &R, mut f: Vec<R::Value>, mut g: Vec<R::Value>) -> Vec<R::Value> {
        if f.is_empty() || g.is_empty() {
            return Vec::new();
        }
        let len = f.len() + g.len() - 1;
        let n = len.next_power_of_two();
        f.resize_with(n, || ring.zero());
        g.resize_with(n, || ring.zero());
        let mut h = convolve_by_transform(&Cyclic::new(n), ring, f, g);
        h.truncate(len);
        h
    }
}

impl<R: RootOfUnity> OnlineConvolution<R> for Additive<Canonical<usize>>
where
    R::Value: Clone,
{
    /// The online polynomial product by the divide and conquer of [`cdq`], each block
    /// contributing through [`convolve_poly`].
    ///
    /// # Complexity
    /// - Time: O(n log^2 n), and `n` calls of `step`
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `R` has no primitive root of unity of order at least `2n`.
    fn online_convolve(
        &self,
        ring: &R,
        g: &[R::Value],
        mut step: impl FnMut(usize, R::Value) -> R::Value,
    ) -> Vec<R::Value> {
        let n = g.len();
        let mut state: (Vec<R::Value>, Vec<R::Value>) = (
            (0..n).map(|_| ring.zero()).collect(),
            (0..n).map(|_| ring.zero()).collect(),
        );
        cdq(
            n,
            &mut state,
            |(f, acc), l, m, r| {
                let h = self.convolve(ring, f[l..m].to_vec(), g[1..r - l].to_vec());
                for t in m..r {
                    acc[t] = ring.add(&acc[t], &h[t - l - 1]);
                }
            },
            |(f, acc), t| f[t] = step(t, acc[t].clone()),
        );
        state.0
    }
}
