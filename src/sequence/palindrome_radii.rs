/// The palindrome radii of `s`.
///
/// # Definition
/// For a sequence `s` of length `n`, `r[i]` is the largest `r` such that `s[i - r + 1..i + r]` is
/// a palindrome, that is the radius of the longest palindrome centered at `i`, counting the center.
///
/// # Complexity
/// - Time: O(n)
/// - Space: O(n)
pub fn palindrome_radii<T: PartialEq>(s: &[T]) -> Vec<usize> {
    let n = s.len();
    let mut r = vec![1; n];
    let (mut c, mut k) = (0, 1);
    for i in 1..n {
        let mut x = if i < k { r[2 * c - i].min(k - i) } else { 0 };
        while x <= i && i + x < n && s[i - x] == s[i + x] {
            x += 1;
        }
        r[i] = x;
        if i + x > k {
            (c, k) = (i, i + x);
        }
    }
    r
}
