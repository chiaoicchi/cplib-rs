use crate::algebra::Monoid;
use crate::range::to_half_open;

/// A sequence `a` of a monoid, with point updates and the folds of its intervals.
///
/// # Definition
/// `n` is the length of `a`. `fold[l, r) = a[l] a[l + 1] ... a[r - 1]`, the product under `op` in
/// this order, which is `id()` for `l = r`.
///
/// # Invariants
/// - `value[n + i] = a[i]` for `i` in `[0, n)`, and `value[i] = op(value[2i], value[2i + 1])` for
///   `i` in `[1, n)`. `value[0]` is unused.
///
/// # Complexity
/// - Space: O(n)
pub struct SegmentTree<M: Monoid> {
    monoid: M,
    value: Box<[M::Value]>,
}

impl<M: Monoid> SegmentTree<M> {
    /// The sequence of `n` copies of `id()`.
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

    /// The segment tree of the sequence `v`.
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

    /// Sets `a[i]` to `x`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i >= n`.
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

    /// Sets `a[i]` to `op(a[i], x)`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i >= n`.
    pub fn op_assign(&mut self, i: usize, x: &M::Value) {
        assert!(
            i < self.len(),
            "index out of bounds: i={i}, len={}",
            self.len()
        );
        self.set(i, self.monoid.op(self.get(i).unwrap(), x));
    }

    /// The element `a[i]`, or `None` if `i >= n`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn get(&self, i: usize) -> Option<&M::Value> {
        self.value.get(self.len().checked_add(i)?)
    }

    /// The fold `fold[l, r)` of `range = [l, r)`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `l > r` or `r > n`.
    pub fn fold(&self, range: impl std::ops::RangeBounds<usize>) -> M::Value {
        let (mut l, mut r) = to_half_open(self.len(), range);
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

    /// An `r` in `[l, n]` at which `pred(fold[l, r))` turns from `true` to `false`.
    ///
    /// # Definition
    /// `pred(fold[l, r))` is `true`, and `r = n` or `pred(fold[l, r + 1))` is `false`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `l > n` or `pred(id())` is `false`.
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

    /// An `l` in `[0, r]` at which `pred(fold[l, r))` turns from `true` to `false` as `l`
    /// decreases.
    ///
    /// # Definition
    /// `pred(fold[l, r))` is `true`, and `l = 0` or `pred(fold[l - 1, r))` is `false`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `r > n` or `pred(id())` is `false`.
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

    /// The length `n`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.value.len() >> 1
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

impl<M: Monoid> std::ops::Index<usize> for SegmentTree<M> {
    type Output = M::Value;
    /// # Panics
    /// Panics if `i >= n`.
    fn index(&self, i: usize) -> &M::Value {
        &self.value[self.len() + i]
    }
}
