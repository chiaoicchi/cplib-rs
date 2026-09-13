/// Returns, for each position of `s`, the previous position holding the same element.
///
/// # Definition
/// For a sequence `s` of length `n`, `p[i]` is `max{j < i: s_j = s_i}` if it exists, and `None`
/// otherwise, for `i` in `[0, n)`. Elements are assumed to be already indexed by `[0, σ)`.
///
/// # Complexity
/// - Time: O(n + σ)
/// - Space: O(n + σ)
///
/// # Panics
/// Panics if `s_i >= sigma`.
pub fn previous_occurrence(s: &[usize], sigma: usize) -> Vec<Option<usize>> {
    let mut p = Vec::with_capacity(s.len());
    let mut last = vec![None; sigma];
    for (i, &c) in s.iter().enumerate() {
        assert!(c < sigma, "element out of range: c={c}, sigma={sigma}");
        p.push(last[c]);
        last[c] = Some(i);
    }
    p
}

/// Returns, for each position of `s`, the next position holding the same element.
///
/// # Definition
/// For a sequence `s` of length `n`, `q[i]` is `min{j > i: s_j = s_i}` if it exists, and `None`
/// otherwise, for `i` in `[0, n)`. Elements are assumed to be already indexed by `[0, σ)`.
///
/// # Complexity
/// - Time: O(n + σ)
/// - Space: O(n + σ)
///
/// # Panics
/// Panics if `s_i >= sigma`.
pub fn next_occurrence(s: &[usize], sigma: usize) -> Vec<Option<usize>> {
    let mut q = Vec::with_capacity(s.len());
    let mut last = vec![None; sigma];
    for (i, &c) in s.iter().enumerate().rev() {
        assert!(c < sigma, "element out of range: c={c}, sigma={sigma}");
        q.push(last[c]);
        last[c] = Some(i);
    }
    q.reverse();
    q
}
