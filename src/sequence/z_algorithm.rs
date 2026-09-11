/// Returns the Z-array of `s`.
///
/// # Definition
/// For a sequence `s` of length `n`, `z[i]` is the length of the longest common prefix of
/// `s` and `s[i..]` for `i` in `[0, n)`.
///
/// # Complexity
/// - Time: O(n)
/// - Space: O(n)
pub fn z_algorithm<T: PartialEq>(s: &[T]) -> Vec<usize> {
    let n = s.len();
    let mut z = vec![n; n];
    let (mut l, mut r) = (0, 1);
    for i in 1..n {
        let mut x = if i < r { z[i - l].min(r - i) } else { 0 };
        while i + x < n && s[x] == s[i + x] {
            x += 1;
        }
        z[i] = x;
        if i + x > r {
            (l, r) = (i, i + x);
        }
    }
    z
}
