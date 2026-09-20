use crate::algebra::Monoid;
use crate::range::to_half_open;

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

    /// Returns an `r` in `[l, n]` at which `pred` switches from `true` to `false`.
    ///
    /// # Definition
    /// Returns some `r` in `[l, n]` such that
    /// - `pred(fold(l..r))` is `true`
    /// - `r = n` or `pred(fold(l..r + 1))` is `false`
    ///
    /// If `pred(fold(l..r))` is monotone in `r` (once `false`, it stays `false`), this is the
    /// maximum `r` such that `pred(fold(l..r))` is `true`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `l > self.len()`.
    /// Panics if `pred(id())` is false.
    pub fn max_right(&self, l: usize, mut pred: impl FnMut(&M::Value) -> bool) -> usize {
        let n = self.len();
        assert!(l <= n, "index out of bounds: l={l}, len={n}");
        let mut acc = self.monoid.id();
        assert!(pred(&acc), "pred(id()) must be true");
        let mut l = l + n;
        let r = n << 1;
        while l < r {
            let k = l.trailing_zeros().min((r - l).ilog2());
            let mut i = l >> k;
            let next = self.monoid.op(&acc, &self.value[i]);
            if pred(&next) {
                acc = next;
                l += 1 << k;
                continue;
            }
            while i < n {
                i <<= 1;
                let next = self.monoid.op(&acc, &self.value[i]);
                if pred(&next) {
                    acc = next;
                    i |= 1;
                }
            }
            return i - n;
        }
        n
    }

    /// Returns an `l` in `[0, r]` at which `pred` switches from `true` to `false`.
    ///
    /// # Definition
    /// Returns some `l` in `[0, r]` such that
    /// - `pred(fold(l..r))` is `true`
    /// - `l = 0` or `pred(fold(l - 1..r))` is `false`.
    ///
    /// If `pred(fold(l..r))` is monotone in `l` (once `false` as `l` decreases, it stays `false`),
    /// this is the minimum `l` such that `pred(fold(l..r))` is `true`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `r > self.len()`.
    /// Panics if `pred(id())` is `false`.
    pub fn min_left(&self, r: usize, mut pred: impl FnMut(&M::Value) -> bool) -> usize {
        let n = self.len();
        assert!(r <= n, "index out of bounds: r={r}, len={n}");
        let mut acc = self.monoid.id();
        assert!(pred(&acc), "pred(id()) must be true");
        let l = n;
        let mut r = r + n;
        while l < r {
            let k = r.trailing_zeros().min((r - l).ilog2());
            let mut i = (r >> k) - 1;
            let next = self.monoid.op(&self.value[i], &acc);
            if pred(&next) {
                acc = next;
                r -= 1 << k;
                continue;
            }
            while i < n {
                i = (i << 1) | 1;
                let next = self.monoid.op(&self.value[i], &acc);
                if pred(&next) {
                    acc = next;
                    i -= 1;
                }
            }
            return i + 1 - n;
        }
        0
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
