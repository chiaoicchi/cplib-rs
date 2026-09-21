use crate::algebra::RootOfUnity;
use crate::divide_and_conquer::cdq;
use crate::poly::poly_convolve;

/// The solution `f` of `f(i) = step(i, Σ_{j=1,...,i} g(j) f(i - j))`, of precision `n`.
///
/// # Definition
/// `f` is the series with `f(i) = step(i, Σ_{j=1,...,i} g(j) f(i - j))` for `0 <= i < n`, the sum
/// being the coefficient of `x^i` in `f g` without the term `g(0) f(i)`. `step` is called for
/// `i = 0, 1, ..., n - 1` in this order.
///
/// # Complexity
/// - Time: O(n log^2 n), and `n` calls of `step`
/// - Space: O(n)
///
/// # Panics
/// Panics if [`poly_convolve`] panics on operands of lengths `[n/2]` and `n - 1`, that is if `R`
/// has no primitive `N`-th root of unity for the least power of two `N >= [n/2] + n - 2`.
pub fn semi_relaxed<R: RootOfUnity<Value: PartialEq + Clone> + Default>(
    ring: &R,
    g: &[R::Value],
    n: usize,
    mut step: impl FnMut(usize, R::Value) -> R::Value,
) -> Vec<R::Value> {
    let mut state: (Vec<R::Value>, Vec<R::Value>) = (
        (0..n).map(|_| ring.zero()).collect(),
        (0..n).map(|_| ring.zero()).collect(),
    );
    cdq(
        n,
        &mut state,
        |(f, acc), l, m, r| {
            let h = poly_convolve(ring, f[l..m].to_vec(), g[1..r - l].to_vec());
            for t in m..r {
                acc[t] = ring.add(&acc[t], &h[t - l - 1]);
            }
        },
        |(f, acc), t| f[t] = step(t, acc[t].clone()),
    );
    state.0
}
