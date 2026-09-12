use crate::divide_and_conquer::cdq;

/// The row minima of a monotone matrix.
///
/// # Definition
/// For an `n x m` matrix `A` given by `f(i, j) = A[i][j]`, returns `j*(i)` for every row `i`,
/// where `j*(i)` is the least `j` attaining `min_j A[i][j]`. `A` is monotone if `j*` is
/// non-decreasing in `i`.
///
/// # Contract
/// `f` is monotone. Total monotonicity suffices: `A[i][j] > A[i][j']` implies
/// `A[i'][j] > A[i'][j']` for `i < i'`, `j < j'`; so does the Monge property
/// `A[i][j] + A[i'][j'] <= A[i][j'] + A[i'][j]`. Unlike monotonicity, these two are inherited by
/// submatrices.
///
/// # Complexity
/// - Time: O(n + m log n)
/// - Space: O(log n)
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

/// `dp[0] = zero` and `dp[i] = min_{j < i} (dp[j] + cost(j, i))` for a Monge `cost`.
///
/// # Definition
/// The matrix `A[i][j] = dp[j] + cost(j, i)`, `j < i` is Monge whenever `cost` is, since the
/// column term `dp[j]` cancels in the quadrangle inequality; hence it is totally monotone, and so
/// is every submatrix.
///
/// # Contract
/// `cost` is Monge: `cost(j, i) + cost(j', i') <= cost(j, i') + cost(j', i)` for
/// `j < j' < i < i'`. Total monotonicity of `cost` alone does not suffice, as the column term
/// `dp[j]` need not preserve it.
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
            |l: usize, m: usize, r: usize| {
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
            |_| {},
        );
    }
    dp.into_iter().map(Option::unwrap).collect()
}
