use crate::algebra::{Action, Monoid};
use crate::range::to_half_open;

/// A sequence `a` of `T`, with the action of a monoid applied on its intervals.
///
/// # Definition
/// `n` is the length of `a`. Applying `f` to `a[i]` replaces it by `act(f, a[i])`.
///
/// # Contract
/// `act` is an action of the monoid on `T`: `act(id(), x) = x` and
/// `act(op(f, g), x) = act(g, act(f, x))` for all `f`, `g` and `x`.
///
/// # Invariants
/// - `a[i] = act(map[p_d], ... act(map[p_1], act(map[p_0], value[i])) ...)`, where `p_0 = n + i`,
///   `p_{t+1} = p_t / 2` and `p_d = 1`. `map[0]` is unused.
///
/// # Complexity
/// - Space: O(n)
pub struct DualSegmentTree<T, M: Monoid, F: Action<T, M::Value>> {
    monoid: M,
    action: F,
    value: Box<[T]>,
    map: Box<[M::Value]>,
}

impl<T, M: Monoid, F: Action<T, M::Value>> DualSegmentTree<T, M, F> {
    /// The dual segment tree of the sequence `v`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    pub fn from_vec(monoid: M, action: F, v: Vec<T>) -> Self {
        let n = v.len();
        Self {
            map: (0..n << 1).map(|_| monoid.id()).collect(),
            monoid,
            action,
            value: v.into(),
        }
    }

    /// Applies `f` to `a[i]`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i >= n`.
    pub fn apply(&mut self, i: usize, f: &M::Value) {
        assert!(
            i < self.len(),
            "index out of bounds: i={i}, len={}",
            self.len()
        );
        self.range_apply(i..=i, f);
    }

    /// Applies `f` to `a[i]` for every `i` in `range = [l, r)`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `l > r` or `r > n`.
    pub fn range_apply(&mut self, range: impl std::ops::RangeBounds<usize>, f: &M::Value) {
        let (mut l, mut r) = to_half_open(self.len(), range);
        l += self.len();
        r += self.len();
        if l == r {
            return;
        }
        self.propagate(l);
        if l != r - 1 {
            self.propagate(r - 1);
        }
        l >>= l.trailing_zeros();
        r >>= r.trailing_zeros();

        loop {
            if l >= r {
                self.map[l] = self.monoid.op(&self.map[l], f);
                l += 1;
                l >>= l.trailing_zeros();
            } else {
                r -= 1;
                self.map[r] = self.monoid.op(&self.map[r], f);
                r >>= r.trailing_zeros();
            }
            if l == r {
                break;
            }
        }
    }

    /// The length `n`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.value.len()
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
    /// the nodes on the path other than `i` hold `id()`.
    fn propagate(&mut self, i: usize) {
        for t in (1..usize::BITS - i.leading_zeros()).rev() {
            let k = i >> t;
            let id = self.monoid.id();
            let f = std::mem::replace(&mut self.map[k], id);
            self.map[k << 1] = self.monoid.op(&self.map[k << 1], &f);
            self.map[(k << 1) | 1] = self.monoid.op(&self.map[(k << 1) | 1], &f);
        }
    }
}

impl<T: Clone, M: Monoid, F: Action<T, M::Value>> DualSegmentTree<T, M, F> {
    /// The element `a[i]`, or `None` if `i >= n`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    pub fn get(&self, mut i: usize) -> Option<T> {
        let mut x = self.value.get(i)?.clone();
        i += self.len();
        while 0 < i {
            x = self.action.act(&self.map[i], &x);
            i >>= 1;
        }
        Some(x)
    }
}
