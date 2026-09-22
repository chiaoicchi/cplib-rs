/// The previous position holding the same element, for each position of `s`.
///
/// # Definition
/// For a sequence `s` of length `n`, `p[i]` is `max{j < i: s[j] = s[i]}` if it exists, and `None`
/// otherwise, for `i` in `[0, n)`.
///
/// # Complexity
/// - Time: O(n + sigma)
/// - Space: O(n + sigma)
///
/// # Panics
/// Panics if `s[i] >= sigma` for some `i`.
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

/// The next position holding the same element, for each position of `s`.
///
/// # Definition
/// For a sequence `s` of length `n`, `q[i]` is `min{j > i: s[j] = s[i]}` if it exists, and `None`
/// otherwise, for `i` in `[0, n)`.
///
/// # Complexity
/// - Time: O(n + sigma)
/// - Space: O(n + sigma)
///
/// # Panics
/// Panics if `s_i >= sigma` for some `i`.
pub fn next_occurrence(s: &[usize], sigma: usize) -> Vec<Option<usize>> {
    let mut q = vec![None; s.len()];
    let mut last = vec![None; sigma];
    for (i, &c) in s.iter().enumerate().rev() {
        assert!(c < sigma, "element out of range: c={c}, sigma={sigma}");
        q[i] = last[c];
        last[c] = Some(i);
    }
    q
}
