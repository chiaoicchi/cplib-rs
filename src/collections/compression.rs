/// A coordinate compression of a finite subset `X` of a totally ordered set `T`.
///
/// # Definition
/// `k = |X|`, and `φ: X -> [0, k)` is the order isomorphism, `φ(x) = #{y in X: y < x}`.
///
/// # Invariants
/// - `value` is `X` in increasing order, so that `value[φ(x)] = x` for `x` in `X`.
///
/// # Complexity
/// - Space: O(k)
pub struct Compression<T> {
    value: Box<[T]>,
}

impl<T: Ord> Compression<T> {
    /// The compression of the set `X` of the elements of `v`.
    ///
    /// # Complexity
    /// - Time: O(n log n), where `n = v.len()`
    /// - Space: O(n)
    pub fn from_vec(mut v: Vec<T>) -> Self {
        v.sort_unstable();
        v.dedup();
        Self { value: v.into() }
    }

    /// `φ(x)`, or `None` if `x` is not in `X`.
    ///
    /// # Complexity
    /// - Time: O(log k)
    /// - Space: O(1)
    pub fn compress(&self, x: &T) -> Option<usize> {
        self.value.binary_search(x).ok()
    }

    /// `#{y in X: y < x}`.
    ///
    /// # Complexity
    /// - Time: O(log k)
    /// - Space: O(1)
    pub fn lower_bound(&self, x: &T) -> usize {
        self.value.partition_point(|y| y < x)
    }

    /// The size `k`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.value.len()
    }

    /// Whether `k = 0`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        self.value.is_empty()
    }
}

impl<T> std::ops::Index<usize> for Compression<T> {
    type Output = T;
    /// The element `φ^{-1}(i)`.
    ///
    /// # Panics
    /// Panics if `i >= k`.
    fn index(&self, i: usize) -> &T {
        assert!(
            i < self.value.len(),
            "index out of bounds: i={i}, len={}",
            self.value.len()
        );
        &self.value[i]
    }
}
