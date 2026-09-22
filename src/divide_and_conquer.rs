/// Calls `contribute` and `finish` over `[0, n)` in the order of the divide and conquer of CDQ.
///
/// # Definition
/// `finish(i)` is called for `i = 0, 1, ..., n - 1` in this order. Interleaved with these,
/// `contribute(l, m, r)` is called with `l < m < r <= n`, so that every pair `i < j` in `[0, n)`
/// lies in exactly one block `[l, m) x [m, r)`, and each call comes after `finish(i)` for every `i`
/// in `[l, m)` and before `finish(j)` for every `j` in `[m, r)`. For `n = 0` nothing is called.
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
