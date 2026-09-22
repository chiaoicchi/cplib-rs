use crate::algebra::{Action, Monoid};
use crate::range::to_half_open;

/// A sequence `a` of a monoid `M`, with the action of a monoid `N` applied on its intervals and
/// the folds of its intervals.
///
/// # Definition
/// `n` is the length of `a`. `fold[l, r) = a[l] a[l + 1] ... a[r - 1]`, the product under `op` of
/// `M` in this order, which is `id()` for `l = r`. Applying `f` to `a[i]` replaces it by
/// `act(f, a[i])`.
///
/// # Contract
/// `act` is an action of `N` on `M` that distributes over `op` of `M`:
/// - `act(map_monoid.id(), x) = x`
/// - `act(map_monoid.op(f, g), x) = act(g, act(f, x))`
/// - `act(f, value_monoid.op(x, y)) = value_monoid.op(act(f, x), act(f, y))`.
///
/// for all `f`, `g`, `x` and `y`.
///
/// # Invariants
/// - `value[k] = act(map[k], op(value[2k, value[2k + 1]))` for `k` in `[1, n)`.
/// - `a[i] = act(map[p_d], ... act(map[p_1], value[p_0]) ...)`, where `p_0 = n + i`,
///   `p_{t+1} = p_t / 2` and `p_d = 1`.
/// - `value[0]` and `map[0]` are unused.
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
    /// The sequence of `n` copies of `id()`.
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

    /// The lazy segment tree of the sequence `v`.
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
        self.propagate(i);
        self.value[i] = x;
        self.pull(i);
    }

    /// Sets `a[i]` to `op(a[i], x)`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i >= n`.
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

    /// Applies `f` to `a[i]`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i >= n`.
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

    /// Applies `f` to `a[i]` for every `i` in `range = [l, r)`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `l > r` or `r > n`.
    pub fn range_apply(&mut self, range: impl std::ops::RangeBounds<usize>, f: &N::Value) {
        let (mut l, mut r) = to_half_open(self.len(), range);
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

    /// The element `a[i]`, or `None` if `i >= n`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    pub fn get(&self, i: usize) -> Option<M::Value> {
        (i < self.len()).then(|| self.fold(i..=i))
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

    /// The length of `n`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.map.len()
    }

    /// Whether `n = 0`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Pushes the maps on the path from the root to the node `i` down to their children, so that
    /// the proper ancestors of `i` hold `id()`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    fn propagate(&mut self, i: usize) {
        for t in (1..usize::BITS - i.leading_zeros()).rev() {
            self.push(i >> t);
        }
    }

    /// Pushes the map of the internal node `k` down to its two children, so that `k` holds `id()`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
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

    /// Recomputes the values of the proper ancestors of the node `i` from their children.
    ///
    /// # Contract
    /// The proper ancestors of `i` hold `id()` in `map`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    fn pull(&mut self, mut i: usize) {
        while 1 < i {
            i >>= 1;
            self.value[i] = self
                .value_monoid
                .op(&self.value[i << 1], &self.value[(i << 1) | 1]);
        }
    }
}
