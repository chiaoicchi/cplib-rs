/// The subsequence automaton of `s`.
///
/// # Definition
/// For a sequence `s` of length `n` over `[0, sigma)`, the subsequence automaton is the deterministic
/// finite automaton with states `0, ..., n`, initial state `0`, all states accepting, and transition
/// `δ(i, c) = min{j >= i: s[j] = c} + 1`, undefined if there is no such `j`.
///
/// # Invariants
/// - `next[i * sigma + c]` is `δ(i, c) - 1`, that is the position of the next `c` at or after `i`,
///   and `u32::MAX` if `δ(i, c)` is undefined, for `i` in `[0, n]` and `c` in `[0, sigma)`.
/// - `len = n`.
///
/// # Complexity
/// - Space: O(n sigma)
pub struct SubsequenceAutomaton {
    next: Box<[u32]>,
    len: usize,
    sigma: usize,
}

impl SubsequenceAutomaton {
    /// The subsequence automaton of `s`.
    ///
    /// # Complexity
    /// - Time: O(n sigma)
    /// - Space: O(n sigma)
    ///
    /// # Panics
    /// Panics if `s[i] >= sigma` for some `i` or `n >= 2^32 - 1`.
    pub fn from_slice(s: &[usize], sigma: usize) -> Self {
        let n = s.len();
        assert!(n < u32::MAX as usize, "n must be less than 2^32 - 1: n={n}");
        let mut next = vec![u32::MAX; (n + 1) * sigma];
        for i in (0..n).rev() {
            let c = s[i];
            assert!(c < sigma, "element out of range: c={c}, sigma={sigma}");
            next.copy_within((i + 1) * sigma..(i + 2) * sigma, i * sigma);
            next[i * sigma + c] = i as u32;
        }
        Self {
            next: next.into_boxed_slice(),
            len: n,
            sigma,
        }
    }

    /// The position `min{j >= i: s[j] = c}` of the next `c` at or after `i`, or `None` if there is
    /// none.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i > n` or `c >= sigma`.
    pub fn next(&self, i: usize, c: usize) -> Option<usize> {
        let n = self.len();
        assert!(i <= n, "state out of range: i={i}, n={n}");
        assert!(
            c < self.sigma,
            "element out of range: c={c}, sigma={}",
            self.sigma
        );
        let j = self.next[i * self.sigma + c];
        (j != u32::MAX).then_some(j as usize)
    }

    /// Whether `t` is a subsequence of `s`.
    ///
    /// # Complexity
    /// - Time: O(|t|)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `t[i] >= sigma` for some `i`.
    pub fn accepts(&self, t: &[usize]) -> bool {
        if let Some(&c) = t.iter().find(|&&c| c >= self.sigma) {
            panic!("element out of range: c={c}, sigma={}", self.sigma);
        }
        let mut i = 0;
        for &c in t {
            match self.next(i, c) {
                Some(j) => i = j + 1,
                None => return false,
            }
        }
        true
    }

    /// The length `n` of `s`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether `n = 0`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The alphabet size `sigma`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn sigma(&self) -> usize {
        self.sigma
    }
}
