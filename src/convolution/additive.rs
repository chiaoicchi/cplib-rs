use crate::algebra::RootOfUnity;
use crate::algebra::additive::Additive;
use crate::algebra::canonical::Canonical;
use crate::convolution::cyclic::convolve_poly;
use crate::convolution::{Convolution, OnlineConvolution};
use crate::divide_and_conquer::cdq;

impl<R: RootOfUnity> Convolution<R> for Additive<Canonical<usize>>
where
    R::Value: Clone,
{
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
        convolve_poly(ring, f.to_vec(), g.to_vec())
    }
}

impl<R: RootOfUnity> OnlineConvolution<R> for Additive<Canonical<usize>>
where
    R::Value: Clone,
{
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
                let h = convolve_poly(ring, f[l..m].to_vec(), g[1..r - l].to_vec());
                for t in m..r {
                    acc[t] = ring.add(&acc[t], &h[t - l - 1]);
                }
            },
            |(f, acc), t| f[t] = step(t, acc[t].clone()),
        );
        state.0
    }
}
