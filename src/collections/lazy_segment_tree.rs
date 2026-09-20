use crate::algebra::{Action, Monoid};
use crate::range::to_half_open;

/// A lazy segment tree data structure.
///
/// # Contract
/// `act` must be an action of the monoid `(N::Value, map_id, map_op)` on the monoid
/// `(M::Value, id, op)`, i.e.
/// `act(map_id(), x) == x` for all `x`.
/// `act(map_op(f, g), x) == act(g, act(f, x))` for all `f`, `g`, `x`.
/// `act(f, op(x, y)) == op(act(f, x), act(f, y))` for all `f`, `x`, `y`.
///
/// # Invariants
/// The values are stored in the internal array `value` as a 1-indexed binary tree in
/// `value[1..2n)`, and the maps in `map[1..n)`; `value[0]` and `map[0]` are unused.
/// - `value[i] = op(value[2i], value[2i + 1])` for all `i` in `[1, n)`, where the map has already
///   been applied to `value[i]` but not to its descendants.
/// - `map[i]` is a map that has not yet been applied to the descendants of `i`.
/// - The map at a node is newer than the maps at its descendants.
///
/// # Complexity
/// - Space: O(n)
pub struct LazySegmentTree<M: Monoid, N: Monoid, F: Action<M::Value, N::Value>> {
    value_monoid: M,
    map_monoid: N,
    action: F,
    value: Box<[M::Value]>,
    map: Box<[N::Value]>,
}

impl<M: Monoid, N: Monoid, F: Action<M::Value, N::Value>> LazySegmentTree<M, N, F> {
    /// Constructs a lazy segment tree with `n` elements, all initialized to `value_monoid.id()`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    pub fn new(value_monoid: M, map_monoid: N, action: F, n: usize) -> Self {
        Self {
            value: (0..n << 1).map(|_| value_monoid.id()).collect(),
            map: (0..n).map(|_| map_monoid.id()).collect(),
            value_monoid,
            map_monoid,
            action,
        }
    }

    /// Constructs a lazy segment tree from a vector `v` of elements.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    pub fn from_vec(value_monoid: M, map_monoid: N, action: F, v: Vec<M::Value>) -> Self {
        let n = v.len();
        let mut value: Vec<M::Value> = (0..n).map(|_| value_monoid.id()).chain(v).collect();
        for i in (1..n).rev() {
            value[i] = value_monoid.op(&value[i << 1], &value[(i << 1) | 1]);
        }
        Self {
            map: (0..n).map(|_| map_monoid.id()).collect(),
            value_monoid,
            map_monoid,
            action,
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
        self.propagate(i);
        self.value[i] = x;
        self.pull(i);
    }

    /// Sets the element at index `i` to `op(a[i], x)`, where `a[i]` is the current element.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i` is out of bounds.
    pub fn op_assign(&mut self, mut i: usize, x: &M::Value) {
        assert!(
            i < self.len(),
            "index out of bounds: i={i}, len={}",
            self.len()
        );
        i += self.len();
        self.propagate(i);
        self.value[i] = self.value_monoid.op(&self.value[i], x);
        self.pull(i);
    }

    /// Sets the element at index `i` to `act(f, a[i])`, where `a[i]` is the current element.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i` is out of bounds.
    pub fn apply(&mut self, mut i: usize, f: &N::Value) {
        assert!(
            i < self.len(),
            "index out of bounds: i={i}, len={}",
            self.len()
        );
        i += self.len();
        self.propagate(i);
        self.value[i] = self.action.act(f, &self.value[i]);
        self.pull(i);
    }

    /// Applies the map `f` to the elements in `range`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if the start of `range` is greater than the end, or the end is greater than
    /// `self.len()`.
    pub fn range_apply(&mut self, range: impl std::ops::RangeBounds<usize>, f: &N::Value) {
        let (mut l, mut r) = to_half_open(self.len(), range);
        assert!(
            l <= r,
            "left bound must be less than or equal to right bound: l={l}, r={r}"
        );
        assert!(r <= self.len(), "range out of bounds: range=[{l}, {r})");
        if l == r {
            return;
        }
        l += self.len();
        r += self.len();
        l >>= l.trailing_zeros();
        r >>= r.trailing_zeros();
        self.propagate(l);
        if l != r - 1 {
            self.propagate(r - 1);
        }
        {
            let (mut l, mut r) = (l, r);
            loop {
                if l >= r {
                    self.value[l] = self.action.act(f, &self.value[l]);
                    if l < self.len() {
                        self.map[l] = self.map_monoid.op(&self.map[l], f);
                    }
                    l += 1;
                    l >>= l.trailing_zeros();
                } else {
                    r -= 1;
                    self.value[r] = self.action.act(f, &self.value[r]);
                    if r < self.len() {
                        self.map[r] = self.map_monoid.op(&self.map[r], f);
                    }
                    r >>= r.trailing_zeros();
                }
                if l == r {
                    break;
                }
            }
        }
        self.pull(l);
        self.pull(r - 1);
    }

    /// Returns a reference to the element at index `i`, or `None` if `i` is out of bounds.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    pub fn get(&mut self, mut i: usize) -> Option<&M::Value> {
        if self.len() <= i {
            return None;
        }
        i += self.len();
        self.propagate(i);
        self.value.get(i)
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
    /// Panics if the start of `range` is greater than the end, or the end is greater than
    /// `self.len()`.
    pub fn fold(&self, range: impl std::ops::RangeBounds<usize>) -> M::Value {
        let (mut l, mut r) = to_half_open(self.len(), range);
        assert!(
            l <= r,
            "left bound must be less than or equal to right bound: l={l}, r={r}"
        );
        assert!(r <= self.len(), "range out of bounds: range=[{l}, {r})");
        if l == r {
            return self.value_monoid.id();
        }
        l += self.len();
        r += self.len();
        l >>= l.trailing_zeros();
        r >>= r.trailing_zeros();

        let mut left = self.value_monoid.id();
        let mut right = self.value_monoid.id();
        loop {
            if l >= r {
                let mut i = l >> 1;
                left = self.value_monoid.op(&left, &self.value[l]);
                l += 1;
                l >>= l.trailing_zeros();
                while i > l >> 1 {
                    left = self.action.act(&self.map[i], &left);
                    i >>= 1;
                }
            } else {
                let mut i = r >> 1;
                r -= 1;
                right = self.value_monoid.op(&self.value[r], &right);
                r >>= r.trailing_zeros();
                while i > r >> 1 {
                    right = self.action.act(&self.map[i], &right);
                    i >>= 1;
                }
            }
            if l == r {
                break;
            }
        }
        let mut x = self.value_monoid.op(&left, &right);
        let mut i = l >> 1;
        while 0 < i {
            x = self.action.act(&self.map[i], &x);
            i >>= 1;
        }
        x
    }

    /// Returns an `r` in `[l, n]` at which `pred` switches from `true` to `false`.
    ///
    /// # Definition
    /// Returns some `r` in `[l, n]` such that
    /// - `pred(fold(l..r))` is `true`
    /// - `r = n` or `pred(fold(l..r + 1))` is `false`
    ///
    /// If `pred(fold(l..r))` is monotone in `r` (one `false`, it stays `false`), this is the
    /// maximum `r` such that `pred(fold(l..r))` is `true`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `l > self.len()`.
    /// Panics if `pred(id())` is `false`.
    pub fn max_right(&mut self, l: usize, mut pred: impl FnMut(&M::Value) -> bool) -> usize {
        let n = self.len();
        assert!(l <= n, "index out of bounds: l={l}, len={n}");
        let mut acc = self.value_monoid.id();
        assert!(pred(&acc), "pred(id()) must be true");
        if l == n {
            return n;
        }
        let mut l = l + n;
        let r = n << 1;
        self.propagate(l);
        self.propagate(r - 1);
        while l < r {
            let k = l.trailing_zeros().min((r - l).ilog2());
            let mut i = l >> k;
            let next = self.value_monoid.op(&acc, &self.value[i]);
            if pred(&next) {
                acc = next;
                l += 1 << k;
                continue;
            }
            while i < n {
                self.push(i);
                i <<= 1;
                let next = self.value_monoid.op(&acc, &self.value[i]);
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
    /// - `l = 0` or `pred(fold(l - 1..r))` is `false`
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
    pub fn min_left(&mut self, r: usize, mut pred: impl FnMut(&M::Value) -> bool) -> usize {
        let n = self.len();
        assert!(r <= n, "index out of bounds: r={r}, len={n}");
        let mut acc = self.value_monoid.id();
        assert!(pred(&acc), "pred(id()) must be true");
        if r == 0 {
            return 0;
        }
        let l = n;
        let mut r = r + n;
        self.propagate(l);
        self.propagate(r - 1);
        while l < r {
            let k = r.trailing_zeros().min((r - l).ilog2());
            let mut i = (r >> k) - 1;
            let next = self.value_monoid.op(&self.value[i], &acc);
            if pred(&next) {
                acc = next;
                r -= 1 << k;
                continue;
            }
            while i < n {
                self.push(i);
                i = (i << 1) | 1;
                let next = self.value_monoid.op(&self.value[i], &acc);
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
        self.map.len()
    }

    /// Returns `true` if the lazy segment tree contains no elements.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Pushes the maps on the path from the root to the leaf `i` down to their children,
    /// so that every node on the path except the leaf holds `id()` afterwards.
    fn propagate(&mut self, i: usize) {
        for t in (1..usize::BITS - i.leading_zeros()).rev() {
            self.push(i >> t);
        }
    }

    /// Pushes the map of the internal node `k` down to its two children, so that `k` holds `id()`
    /// afterwards.
    fn push(&mut self, k: usize) {
        let id = self.map_monoid.id();
        let f = std::mem::replace(&mut self.map[k], id);
        self.value[k << 1] = self.action.act(&f, &self.value[k << 1]);
        self.value[(k << 1) | 1] = self.action.act(&f, &self.value[(k << 1) | 1]);
        if k << 1 < self.len() {
            self.map[k << 1] = self.map_monoid.op(&self.map[k << 1], &f);
        }
        if (k << 1) + 1 < self.len() {
            self.map[(k << 1) | 1] = self.map_monoid.op(&self.map[(k << 1) | 1], &f);
        }
    }

    /// Recomputes the values of the ancestors of the leaf `i` from their children.
    fn pull(&mut self, mut i: usize) {
        while 1 < i {
            i >>= 1;
            self.value[i] = self
                .value_monoid
                .op(&self.value[i << 1], &self.value[(i << 1) | 1]);
        }
    }
}
