/// The subsequence automaton of `s`.
///
/// # Definition
/// For a sequence `s` of length `n` over `[0, σ)`, the subsequence automaton is the deterministic
/// finite automaton with states `0, ..., n`, initial state `0`, all states accepting, and transition
/// `δ(i, c) = min{j >= i: s_j = c} + 1` (undefined if no such `j`). It accepts exactly the
/// subsequences of `s`, and running it on `t` embeds `t` into `s` leftmost.
///
/// # Invariants
/// `next[i * sigma + c]` is `δ(i, c) - 1`, i.e. the position of the next `c` at or after `i`,
/// and `u32::MAX` if undefined for `i` in `[0, n]` and `c` in `[0, σ)`.
/// `len` is `n` and `sigma` is `σ`.
///
/// # Complexity
/// - Space: O(nσ)
pub struct SubsequenceAutomaton {
    next: Box<[u32]>,
    len: usize,
    sigma: usize,
}

impl SubsequenceAutomaton {
    /// Constructs the subsequence automaton of `s`.
    ///
    /// # Complexity
    /// - Time: O(nσ)
    /// - Space: O(nσ)
    ///
    /// # Panics
    /// Panics if `s_i >= sigma`.
    /// Panics if `n >= 2^32 - 1`.
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

    /// Returns `min{j >= i: s_j = c}`, and `None` if no such `j` exists.
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

    /// Returns `true` if `t` is a subsequence of `s`.
    ///
    /// # Complexity
    /// - Time: O(|t|)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `t_i >= sigma`.
    pub fn accepts(&self, t: &[usize]) -> bool {
        let mut i = 0;
        for &c in t {
            match self.next(i, c) {
                Some(j) => i = j + 1,
                None => return false,
            }
        }
        true
    }

    /// Returns `n`, the length of `s`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.len
    }

    /// Returns `true` if `s` is empty.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns `σ`, the alphabet size.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn sigma(&self) -> usize {
        self.sigma
    }
}
