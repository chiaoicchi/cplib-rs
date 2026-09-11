/// Returns the border array of `s`.
///
/// # Definition
/// For a sequence `s` of length `n`, `b[i]` is the length of the longest border of `s[0..i)`,
/// i.e. the longest proper prefix of `s[0..i)` that is also its suffix, for `i` in `[0, n]`.
/// The borders of `s[0..i)` are exactly `b[i]`, `b[b[i]]`, ..., down to `0`, and the smallest
/// period of `s[0..i)` is `i - b[i]`.
///
/// # Complexity
/// - Time: O(n)
/// - Space: O(n)
pub fn border_array<T: PartialEq>(s: &[T]) -> Vec<usize> {
    let n = s.len();
    let mut b = vec![0; n + 1];
    for i in 2..=n {
        let mut x = b[i - 1];
        while x > 0 && s[i - 1] != s[x] {
            x = b[x];
        }
        if s[i - 1] == s[x] {
            x += 1;
        }
        b[i] = x;
    }
    b
}
