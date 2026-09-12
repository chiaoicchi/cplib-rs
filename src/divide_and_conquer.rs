/// Runs the divide and conquer of CDQ on `[0, n)`.
///
/// # Definition
/// For `[l, r)` with `r - l >= 2` and `m = l + (r - l) / 2`, `cdq` processes `[l, m)`, calls
/// `contribute(l, m, r)`, then processes `[m, r)`; for `r - l = 1` it calls `finish(l)`. Hence
/// every pair `i < j` in `[0, n)` is covered by exactly one call `contribute(l, m, r)` with
/// `l <= i < m <= j < r`, all such calls for a fixed `j` precede `finish(j)`, and `finish(i)`
/// precedes every call reading `i`.
///
/// # Contract
/// - `contribute(l, m, r)` reads only indices in `[l, m)` and writes only indices in `[m, r)`.
/// - `finish(i)` reads and writes only index `i`.
///
/// # Complexity
/// - Time: O(n) calls of `contribute`, whose widths `r - l` sum to O(n log n), and `n` calls of
///   `finish`
/// - Space: O(log n)
pub fn cdq<S>(
    n: usize,
    state: &mut S,
    mut contribute: impl FnMut(&mut S, usize, usize, usize),
    mut finish: impl FnMut(&mut S, usize),
) {
    fn rec<S>(
        l: usize,
        r: usize,
        state: &mut S,
        contribute: &mut impl FnMut(&mut S, usize, usize, usize),
        finish: &mut impl FnMut(&mut S, usize),
    ) {
        if r - 1 == l {
            finish(state, l);
            return;
        }
        let m = l + (r - l) / 2;
        rec(l, m, state, contribute, finish);
        contribute(state, l, m, r);
        rec(m, r, state, contribute, finish);
    }
    if n > 0 {
        rec(0, n, state, &mut contribute, &mut finish);
    }
}
