//! Monotone and Monge matrices.
//!
//! For a matrix `A`, `j*(i)` is the least `j` attaining `min_j A[i][j]`. `A` is monotone if `j*` is
//! non-decreasing. `A` is totally monotone if `A[i][j] > A[i][j']` implies `A[i'][j] > A[i'][j']`,
//! and Monge if `A[i][j] + A[i'][j'] <= A[i][j'] + A[i'][j]`, both for all `i < i'` and `j < j'` at
//! which the entries are defined.

use crate::divide_and_conquer::cdq;

/// The leftmost minimum of each row of a monotone matrix, as its column.
///
/// # Definition
/// For the `n x m` matrix `A[i][j] = f(i, j)`, returns `j*(0), ..., j*(n - 1)`.
///
/// # Contract
/// `A` is monotone.
///
/// # Complexity
/// - Time: O(n + m log n)
/// - Space: O(n)
///
/// # Panics
/// Panics if `m == 0` and `n > 0`.
pub fn monotone_minima<T: PartialOrd>(
    n: usize,
    m: usize,
    f: impl Fn(usize, usize) -> T,
) -> Vec<usize> {
    assert!(
        m > 0 || n == 0,
        "a nonempty matrix must have a column: n={n}, m={m}"
    );
    fn rec<T: PartialOrd>(
        f: &impl Fn(usize, usize) -> T,
        argmin: &mut [usize],
        (l, r): (usize, usize),
        (cl, cr): (usize, usize),
    ) {
        if l >= r {
            return;
        }
        let o = l + (r - l) / 2;
        let mut best = cl;
        let mut value = f(o, cl);
        for j in cl + 1..=cr {
            let v = f(o, j);
            if v < value {
                (best, value) = (j, v);
            }
        }
        argmin[o] = best;
        rec(f, argmin, (l, o), (cl, best));
        rec(f, argmin, (o + 1, r), (best, cr));
    }
    let mut argmin = vec![0; n];
    if n > 0 {
        rec(&f, &mut argmin, (0, n), (0, m - 1));
    }
    argmin
}

/// The shortest distances from `0` in the DAG on `[0, n)` with the edges `j -> i` for `j < i` of a
/// Monge cost.
///
/// # Definition
/// `dp[0] = zero` and `dp[i] = min_{j<i} (dp[j] + cost(j, i))` for `1 <= i < n`. Returns
/// `dp[0], ..., dp[n - 1]`, which is empty for `n = 0`.
///
/// # Contract
/// The matrix `C[i][j] = cost(j, i)`, defined for `j < i`, is Monge.
///
/// # Complexity
/// - Time: O(n log^2 n)
/// - Space: O(n)
pub fn monge_dp<T: PartialOrd + Clone + std::ops::Add<Output = T>>(
    n: usize,
    zero: T,
    cost: impl Fn(usize, usize) -> T,
) -> Vec<T> {
    let mut dp = vec![None; n];
    if n > 0 {
        dp[0] = Some(zero);
        cdq(
            n,
            &mut dp,
            |dp, l, m, r| {
                let argmin = monotone_minima(r - m, m - l, |i, j| {
                    dp[l + j].clone().unwrap() + cost(l + j, m + i)
                });
                for (i, &j) in argmin.iter().enumerate() {
                    let v = dp[l + j].clone().unwrap() + cost(l + j, m + i);
                    let cur = &mut dp[m + i];
                    if cur.as_ref().is_none_or(|c| v < *c) {
                        *cur = Some(v);
                    }
                }
            },
            |_, _| {},
        );
    }
    dp.into_iter().map(Option::unwrap).collect()
}
