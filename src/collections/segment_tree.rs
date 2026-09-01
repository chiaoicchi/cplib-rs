use crate::algebra::Monoid;

/// A segment tree data structure.
///
/// # Invariants
/// The elements are stored in the internal array `value` as a 1-indexed binary tree in `value[1..2n)`; `value[0]` is unused.
/// - `value[n + i] = a[i]` for all `i` in `[0, n)` (leaves)
/// - `value[i] = op(value[2i], value[2i + 1])` for all `i` in `[1, n)` (internal nodes)
///
/// where `a` denotes the logical sequence of elements.
///
/// # Complexity
/// - Space: O(n)
pub struct SegmentTree<M: Monoid> {
    monoid: M,
    value: Box<[M::Value]>,
}

impl<M: Monoid> SegmentTree<M> {
    /// Constructs a segment tree with `n` elements, all initialized to `id()`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    pub fn new(monoid: M, n: usize) -> Self {
        Self {
            value: (0..n << 1).map(|_| monoid.id()).collect(),
            monoid,
        }
    }

    /// Constructs a segment tree from a vector `v` of elements.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    pub fn from_vec(monoid: M, v: Vec<M::Value>) -> Self {
        let n = v.len();
        let mut value: Vec<M::Value> = (0..n).map(|_| monoid.id()).chain(v).collect();
        for i in (1..n).rev() {
            value[i] = monoid.op(&value[i << 1], &value[(i << 1) | 1]);
        }
        Self {
            monoid,
            value: value.into(),
        }
    }

    /// Sets the element at index `i` to `x`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i` is out of bounds.
    pub fn set(&mut self, mut i: usize, x: M::Value) {
        assert!(
            i < self.len(),
            "index out of bounds: i={i}, len={}",
            self.len()
        );
        i += self.len();
        self.value[i] = x;
        while 1 < i {
            i >>= 1;
            self.value[i] = self
                .monoid
                .op(&self.value[i << 1], &self.value[(i << 1) | 1]);
        }
    }

    /// Sets the element at index `i` to `op(a[i], x)`, where `a[i]` is the current element.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i` is out of bounds.
    pub fn op_assign(&mut self, i: usize, x: &M::Value) {
        assert!(
            i < self.len(),
            "index out of bounds: i={i}, len={}",
            self.len()
        );
        self.set(i, self.monoid.op(self.get(i).unwrap(), x));
    }

    /// Returns a reference to the element at index `i`, or `None` if `i` is out of bounds.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn get(&self, i: usize) -> Option<&M::Value> {
        self.value.get(self.len().checked_add(i)?)
    }

    /// Folds the elements in `range`.
    ///
    /// Returns `op(...op(op(a[l], a[l + 1]), a[l + 2])..., a[r - 1])` where `range` is `l..r`,
    /// or `id()` if `range` is empty.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if the start of `range` is greater than the end, or the end is greater
    /// than `self.len()`.
    pub fn fold(&self, range: impl std::ops::RangeBounds<usize>) -> M::Value {
        let (mut l, mut r) = to_half_open(self.len(), range);
        assert!(
            l <= r,
            "left bound must be less than or equal to right bound: l={l}, r={r}"
        );
        assert!(r <= self.len(), "range out of bounds: range=[{l}, {r})");
        l += self.len();
        r += self.len();
        if l == r {
            return self.monoid.id();
        }
        l >>= l.trailing_zeros();
        r >>= r.trailing_zeros();

        let mut left = self.monoid.id();
        let mut right = self.monoid.id();
        loop {
            if l >= r {
                left = self.monoid.op(&left, &self.value[l]);
                l += 1;
                l >>= l.trailing_zeros();
            } else {
                r -= 1;
                right = self.monoid.op(&self.value[r], &right);
                r >>= r.trailing_zeros();
            }
            if l == r {
                break;
            }
        }
        self.monoid.op(&left, &right)
    }

    /// Returns the number of elements.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.value.len() >> 1
    }

    /// Returns `true` if the segment tree contains no elements.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl<M: Monoid> std::ops::Index<usize> for SegmentTree<M> {
    type Output = M::Value;
    /// Returns a reference to the element at index `i`.
    ///
    /// # Panics
    /// Panics if `i` is out of bounds.
    fn index(&self, i: usize) -> &M::Value {
        &self.value[self.len() + i]
    }
}

/// Converts `range` to a pair `(l, r)` representing the half-open interval `[l, r)`,
/// where an unbounded end defaults to `0` or `max`.
fn to_half_open(max: usize, range: impl std::ops::RangeBounds<usize>) -> (usize, usize) {
    use std::ops::Bound;
    let l = match range.start_bound() {
        Bound::Unbounded => 0,
        Bound::Included(&x) => x,
        Bound::Excluded(&x) => x + 1,
    };
    let r = match range.end_bound() {
        Bound::Unbounded => max,
        Bound::Included(&x) => x + 1,
        Bound::Excluded(&x) => x,
    };
    (l, r)
}
