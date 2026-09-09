/// A coordinate compression of a finite subset `X` of a totally ordered set `T`.
///
/// # Definition
/// Let `X` be a finite subset of `T` with `|X| = k`. There is a unique order isomorphism
/// `φ: X -> [0, k)`, given by `φ(x) = #{y in X: y < x}`.
///
/// # Invariants
/// `value` is strictly increasing and enumerates `X`, so `value[φ(x)] = x` for `x` in `X`.
///
/// # Complexity
/// - Space: O(k)
pub struct Compression<T> {
    value: Box<[T]>,
}

impl<T: Ord> Compression<T> {
    /// Constructs a compression of the set of elements of `v`.
    ///
    /// # Complexity
    /// - Time: O(n log n)
    /// - Space: O(n)
    pub fn from_vec(mut v: Vec<T>) -> Self {
        v.sort_unstable();
        v.dedup();
        Self { value: v.into() }
    }

    /// Returns `φ(x)` if `x` is in `X`, and `None` otherwise.
    ///
    /// # Complexity
    /// - Time: O(log k)
    /// - Space: O(1)
    pub fn compress(&self, x: &T) -> Option<usize> {
        self.value.binary_search(x).ok()
    }

    /// Returns `#{y in X: y < x}`, the extension of `φ` to `T`.
    ///
    /// # Complexity
    /// - Time: O(log k)
    /// - Space: O(1)
    pub fn lower_bound(&self, x: &T) -> usize {
        self.value.partition_point(|y| y < x)
    }

    /// Returns `k = |X|`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.value.len()
    }

    /// Returns `true` if `X` is empty.
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
    /// Returns `φ^{-1}(i)`, the `i`-th smallest element of `X`.
    ///
    /// # Panics
    /// Panics if `i >= self.len()`.
    fn index(&self, i: usize) -> &T {
        assert!(
            i < self.value.len(),
            "index out of bounds: i={i}, len={}",
            self.value.len()
        );
        &self.value[i]
    }
}
