use crate::algebra::RootOfUnity;
use crate::poly::subproduct_tree::SubproductTree;

/// The coefficients of the polynomial `f` of degree less than `n` with `f(x(j)) = y(j)`.
///
/// # Definition
/// For `n = xs.len()`, with `x(j)` and `y(j)` the `j`-th entries of `xs` and `ys`. Returns
/// `f = Σ_j y(j) Π_{i!=j} (t - x(j)) / Π_{i!=j} (x(j) - x(i))`, of length `n`.
///
/// # Contract
/// The points `x(j)` are distinct.
///
/// # Complexity
/// - Time: O(n log^2 n)
/// - Space: O(n log n)
///
/// # Panics
/// Panics if `ys.len() != n`.
/// Panics if `R` has no primitive `N`-th root of unity, where `N` is the least power of two at
/// least `n` if `n <= 32` and at least `2n - 1` otherwise.
pub fn poly_interpolate<R: RootOfUnity + Clone>(
    ring: &R,
    xs: &[R::Value],
    ys: &[R::Value],
) -> Vec<R::Value> {
    SubproductTree::new(ring.clone(), xs).interpolate(ys)
}
