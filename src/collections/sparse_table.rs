use crate::algebra::{Idempotent, Semigroup};
use crate::range::to_half_open;

/// A sparse table, a static structure for range folds of an idempotent semigroup.
///
/// # Definition
/// Let `a` be the sequence of elements of an idempotent semigroup `(S, op)`. For a half-open
/// interval `[l, r)` with `l < r`, `fold[l, r) = a[l] op ... op a[r - 1]`.
///
/// # Invariants
/// `table[k][i] = fold[i, i + 2^k)` for all `k` in `[0, floor(log_2 n)]` and `i` in
/// `[0, n - 2^k]`; in particular `table[0] = a`.
///
/// # Complexity
/// - Space: O(n log n)
pub struct SparseTable<S: Semigroup + Idempotent> {
    semigroup: S,
    table: Box<[Box<[S::Value]>]>,
}

impl<S: Semigroup + Idempotent> SparseTable<S> {
    /// Constructs a sparse table from a vector `v` of elements.
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

    /// Returns the element at index `i`, or `None` if `i` is out of bounds.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn get(&self, i: usize) -> Option<&S::Value> {
        self.table.first()?.get(i)
    }

    /// Returns the fold of the elements in `range`, or `None` if `range` is empty.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `range` is out of bounds or `l > r`.
    pub fn fold(&self, range: impl std::ops::RangeBounds<usize>) -> Option<S::Value> {
        let (l, r) = to_half_open(self.len(), range);
        assert!(
            l <= r,
            "left bound must be less than or equal to right bound: l={l}, r={r}"
        );
        assert!(r <= self.len(), "range out of bounds: range=[{l}, {r})");
        if l == r {
            return None;
        }
        let k = (r - l).ilog2() as usize;
        let row = &self.table[k];
        Some(self.semigroup.op(&row[l], &row[r - (1 << k)]))
    }

    /// Returns the number of elements.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.table.first().map_or(0, |row| row.len())
    }

    /// Returns `true` if the sparse table has no elements.
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
    /// Returns a reference to the element at index `i`.
    ///
    /// # Panics
    /// Panics if `i` is out of bounds.
    fn index(&self, i: usize) -> &S::Value {
        assert!(
            i < self.len(),
            "index out of bounds: i={i}, len={}",
            self.len()
        );
        &self.table.first().unwrap()[i]
    }
}
