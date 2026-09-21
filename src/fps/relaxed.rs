use crate::algebra::RootOfUnity;
use crate::divide_and_conquer::cdq;
use crate::poly::poly_convolve;

/// The series `f` in `R[[x]]/(x^n)` with `f(i) = step(i, Σ_{j=1,...,i} g(j) f(i - j))`, as its `n`
/// coefficients.
///
/// # Definition
/// The sum is the coefficient of `x^i` in `fg` without the term `g(0) f(i)`, and is `0` for
/// `i = 0`. `step` is called for `i = 0, 1, ..., n - 1` in this order.
///
/// # Complexity
/// - Time: O(n log^2 n)
/// - Space: O(n)
///
/// # Panics
/// Panics if `n >= 2` and `g.len() < n`.
/// Panics only if `R` has no primitive `N`-th root of unity, where `N` is the least power of two at
/// least `2n`, apart from the above.
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
