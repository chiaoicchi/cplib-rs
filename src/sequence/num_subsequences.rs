use crate::algebra::One;

/// Returns the number of distinct subsequences of `s`, represented in `N`.
///
/// # Definition
/// A subsequence of `s = (s_0, ..., s_{n-1})` is a sequence `(s_{i_1}, ..., s_{i_k})` with
/// `i_1 < ... < i_k`. Two subsequences are identified if they are equal as sequences.
/// The result is the image in `N` of the integer `#{t: t is a subsequence of s}`. Elements are
/// assumed to be already indexed by `[0, σ)`.
///
/// # Complexity
/// - Time: O(n + σ)
/// - Space: O(σ)
///
/// # Panics
/// Panics if `s_i >= sigma`.
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
