use crate::algebra::min::Min;
use crate::collections::sparse_table::SparseTable;

/// A suffix array with its inverse and the LCP array.
///
/// # Definition
/// Let `s` be a sequence of length `n` over a totally ordered alphabet, and `s[i..]` its suffix
/// starting at `i`. `sa` is the permutation of `[0, n)` such that `s[sa[0]..] < s[sa[1]..] < ...`
/// in lexicographic order, and `isa` is its inverse. `lcp[k]` is the length of the longest common
/// prefix of `s[sa[k]..]` and `s[sa[k + 1]..]` for `k` in `[0, n - 1)`.
///
/// For `i != j` with `isa[i] < isa[j]`, the longest common prefix of `s[i..]` and `s[j..]` has
/// length `min(lcp[isa[i]..isa[j]])`, since the longest common prefix of two suffixes is the
/// minimum over the suffixes between them in lexicographic order.
///
/// # Invariants
/// - `isa[sa[k]] = k` for all `k`.
/// - `lcp.len() = n - 1` if `n > 0` and `0` otherwise.
/// - `rmq`, once built, is a sparse table over `lcp` under `min`.
///
/// # Complexity
/// - Space: O(n), plus O(n log n) once `lcp_of` has been called.
pub struct SuffixArray {
    sa: Box<[usize]>,
    isa: Box<[usize]>,
    lcp: Box<[usize]>,
    rmq: std::cell::OnceCell<SparseTable<Min<usize>>>,
}

impl SuffixArray {
    /// Constructs the suffix array of `s`.
    ///
    /// # Complexity
    /// - Time: O(n log n)
    /// - Space: O(n)
    pub fn from_slice<T: Ord>(s: &[T]) -> Self {
        let n = s.len();
        let mut sa: Vec<usize> = (0..n).collect();
        sa.sort_unstable_by_key(|&i| &s[i]);
        let mut isa = vec![0; n];
        let mut classes = 0;
        for k in 0..n {
            if k == 0 || s[sa[k - 1]] != s[sa[k]] {
                classes += 1;
            }
            isa[sa[k]] = classes;
        }
        let mut tmp = vec![0; n];
        let mut count = vec![0; n + 1];
        let mut k = 1;
        while classes < n && k < n {
            let mut t = 0;
            for i in n - k..n {
                tmp[t] = i;
                t += 1;
            }
            for &j in &sa {
                if j >= k {
                    tmp[t] = j - k;
                    t += 1;
                }
            }
            count.fill(0);
            for &i in &tmp {
                count[isa[i]] += 1;
            }
            let mut sum = 0;
            for c in &mut count {
                sum += *c;
                *c = sum - *c;
            }
            for &i in &tmp {
                sa[count[isa[i]]] = i;
                count[isa[i]] += 1;
            }
            let key = |i: usize| (isa[i], if i + k < n { isa[i + k] } else { 0 });
            classes = 0;
            for j in 0..n {
                if j == 0 || key(sa[j - 1]) != key(sa[j]) {
                    classes += 1;
                }
                tmp[sa[j]] = classes;
            }
            std::mem::swap(&mut isa, &mut tmp);
            k <<= 1;
        }
        for r in &mut isa {
            *r -= 1;
        }

        let mut lcp = vec![0; n.saturating_sub(1)];
        let mut h = 0;
        for i in 0..n {
            if isa[i] == 0 {
                h = 0;
                continue;
            }
            let j = sa[isa[i] - 1];
            while i + h < n && j + h < n && s[i + h] == s[j + h] {
                h += 1;
            }
            lcp[isa[i] - 1] = h;
            h = h.saturating_sub(1);
        }
        Self {
            sa: sa.into(),
            isa: isa.into(),
            lcp: lcp.into(),
            rmq: std::cell::OnceCell::new(),
        }
    }

    /// Returns `sa`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn sa(&self) -> &[usize] {
        &self.sa
    }

    /// Returns `isa`, the inverse of `sa`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn isa(&self) -> &[usize] {
        &self.isa
    }

    /// Returns `lcp`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn lcp(&self) -> &[usize] {
        &self.lcp
    }

    /// Returns the length of the longest common prefix of `s[i..]` and `s[j..]`.
    ///
    /// # Complexity
    /// - Time: O(1), plus O(n log n) on the first call
    /// - Space: O(1), plus O(n log n) on the first call
    ///
    /// # Panics
    /// Panics if `i >= n` or `j >= n`.
    pub fn lcp_of(&self, i: usize, j: usize) -> usize {
        assert!(
            i < self.len(),
            "index out of bounds: i={i}, len={}",
            self.len()
        );
        assert!(
            j < self.len(),
            "index out of bounds: j={j}, len={}",
            self.len()
        );
        if i == j {
            return self.len() - i;
        }
        let (a, b) = (self.isa[i].min(self.isa[j]), self.isa[i].max(self.isa[j]));
        let rmq = self
            .rmq
            .get_or_init(|| SparseTable::from_vec(Min::new(), self.lcp.to_vec()));
        rmq.fold(a..b).unwrap()
    }

    /// Returns the length of the sequence.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.sa.len()
    }

    /// Returns `true` if the sequence is empty.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        self.sa.is_empty()
    }
}
