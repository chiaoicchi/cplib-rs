use crate::algebra::One;

/// The number of distinct subsequences of `s`, including the empty one, in `N`.
///
/// # Definition
/// A subsequence of `s` is a sequence `(s[i_1], ..., s[i_k])` with `k >= 0` and
/// `i_1 < ... < i_k`, and two subsequences are identified if they are equal as sequences.
/// The result is the image in `N` of the number of subsequences of `s`.
///
/// # Complexity
/// - Time: O(n + sigma)
/// - Space: O(sigma)
///
/// # Panics
/// Panics if `s[i] >= sigma` for some `i`.
pub fn num_subsequences<N: Clone + One + std::ops::Add<Output = N> + std::ops::Sub<Output = N>>(
    s: &[usize],
    sigma: usize,
) -> N {
    let mut f = N::one();
    let mut last: Vec<Option<N>> = vec![None; sigma];
    for &c in s {
        assert!(c < sigma, "element out of range: c={c}, sigma={sigma}");
        let g = f.clone() + f.clone();
        let g = match &last[c] {
            Some(x) => g - x.clone(),
            None => g,
        };
        last[c] = Some(std::mem::replace(&mut f, g));
    }
    f
}
