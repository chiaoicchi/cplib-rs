use crate::collections::range_map::RangeMap;

/// A subset of a totally ordered set `K` that is a finite union of half-open intervals,
/// stored as its decomposition into maximal intervals.
///
/// # Definition
/// A finite subset `S` of `K` of the form `[a_1, b_1) or ... or [a_m, b_m)` with
/// `a_1 < b_1 < a_2 < b_2 < ... < a_m < b_m`. The decomposition is unique.
///
/// # Invariants
/// The set is a [`RangeMap<K, ()>`] whose intervals are pairwise non-adjacent, i.e. `b_i < a_j`
/// for consecutive keys. Every mutation restores this by merging.
///
/// # Complexity
/// - Space: O(m), where `m` is the number of intervals
pub struct RangeSet<K>(RangeMap<K, ()>);

impl<K: Ord + Copy> RangeSet<K> {
    /// Constructs an empty set.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn new() -> Self {
        Self(RangeMap::new())
    }

    /// Returns `true` if `p` is in `S`.
    ///
    /// # Complexity
    /// - Time: O(log m)
    /// - Space: O(1)
    pub fn contains(&self, p: K) -> bool {
        self.get(p).is_some()
    }

    /// Returns the maximal interval `[a, b)` of `S` containing `p`, and `None` if `p` is not in
    /// `S`.
    ///
    /// # Complexity
    /// - Time: O(log m)
    /// - Space: O(1)
    pub fn get(&self, p: K) -> Option<(K, K)> {
        let (a, b, _) = self.0.get(p)?;
        Some((a, b))
    }

    /// Adds `[l, r)` to `S`.
    ///
    /// # Complexity
    /// - Time: amortized O(log m)
    /// - Space: amortized O(1)
    ///
    /// # Panics
    /// Panics if `l >= r`.
    pub fn insert(&mut self, l: K, r: K) {
        self.0.assign(l, r, ());
        self.0.merge(l);
    }

    /// Removes `[l, r)` from `S`.
    ///
    /// # Complexity
    /// - Time: amortized O(log m)
    /// - Space: amortized O(1)
    ///
    /// # Panics
    /// Panics if `l >= r`.
    pub fn remove(&mut self, l: K, r: K) {
        self.0.remove(l, r);
    }

    /// Iterates over the maximal intervals intersecting `[l, r)` in increasing order.
    /// Intervals are not clipped.
    ///
    /// # Complexity
    /// - Time: O((k + 1) log m) to exhaust, where `k` is the number of intervals yielded
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `l >= r`.
    pub fn range(&self, l: K, r: K) -> impl DoubleEndedIterator<Item = (K, K)> {
        self.0.range(l, r).map(|(a, b, _)| (a, b))
    }

    /// Iterates over the maximal intervals in increasing order.
    ///
    /// # Complexity
    /// - Time: O(m) to exhaust
    /// - Space: O(1)
    pub fn iter(&self) -> impl DoubleEndedIterator<Item = (K, K)> {
        self.0.iter().map(|(a, b, _)| (a, b))
    }

    /// Returns the number of maximal intervals `m`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns `true` if `S` is empty.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        self.0.len() == 0
    }
}

impl<K: Ord + Copy> Default for RangeSet<K> {
    fn default() -> Self {
        Self::new()
    }
}
