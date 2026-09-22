use crate::algebra::{Idempotent, Semigroup};
use crate::range::to_half_open;

/// A fixed sequence `a` of an idempotent semigroup, with the folds of its intervals.
///
/// # Definition
/// `n` is the length of `a`. For `l < r`, `fold[l, r) = a[l] a[l + 1] ... a[r - 1]`, the product
/// under `op` in this order.
///
/// # Invariants
/// - `table[k][i] = fold[i, i + 2^k)` for `k` in `[0, floor(log2 n)]` and `i` in `[0, n - 2^k]`; in
///   particular `table[0] = a`. `table` is empty for `n = 0`.
///
/// # Complexity
/// - Space: O(n log n)
pub struct SparseTable<S: Semigroup + Idempotent> {
    semigroup: S,
    table: Box<[Box<[S::Value]>]>,
}

impl<S: Semigroup + Idempotent> SparseTable<S> {
    /// The sparse table of the sequence `v`.
    ///
    /// # Complexity
    /// - Time: O(n log n)
    /// - Space: O(n log n)
    pub fn from_vec(semigroup: S, v: Vec<S::Value>) -> Self {
        let n = v.len();
        let mut table: Vec<Box<[S::Value]>> = Vec::new();
        if n > 0 {
            table.push(v.into());
            for k in 1..=n.ilog2() as usize {
                let half = 1 << (k - 1);
                let prev = &table[k - 1];
                let next: Box<[S::Value]> = (0..=n - (1 << k))
                    .map(|i| semigroup.op(&prev[i], &prev[i + half]))
                    .collect();
                table.push(next);
            }
        }
        Self {
            semigroup,
            table: table.into(),
        }
    }

    /// The element `a[i]`, or `None` if `i >= n`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn get(&self, i: usize) -> Option<&S::Value> {
        self.table.first()?.get(i)
    }

    /// The fold `fold[l, r)` of `range = [l, r)`, or `None` if `l = r`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `l > r` or `r > n`.
    pub fn fold(&self, range: impl std::ops::RangeBounds<usize>) -> Option<S::Value> {
        let (l, r) = to_half_open(self.len(), range);
        if l == r {
            return None;
        }
        let k = (r - l).ilog2() as usize;
        let row = &self.table[k];
        Some(self.semigroup.op(&row[l], &row[r - (1 << k)]))
    }

    /// The length `n`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.table.first().map_or(0, |row| row.len())
    }

    /// Whether `n = 0`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl<S: Semigroup + Idempotent> std::ops::Index<usize> for SparseTable<S> {
    type Output = S::Value;
    /// # Panics
    /// Panics if `i >= n`.
    fn index(&self, i: usize) -> &S::Value {
        assert!(
            i < self.len(),
            "index out of bounds: i={i}, len={}",
            self.len()
        );
        &self.table.first().unwrap()[i]
    }
}
