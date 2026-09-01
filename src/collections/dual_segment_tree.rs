use crate::algebra::{Action, Monoid};
use crate::range::to_half_open;

/// A dual segment tree data structure.
///
/// # Contract
/// `act` must be an action of the monoid `(M::Value, id, op)` on `T`, i.e.
/// `act(id(), x) == x` and `act(op(f, g), x) == act(g, act(f, x))` for all `f`, `g`, `x`.
///
/// # Invariants
/// The maps are stored in the internal array `map` as a 1-indexed binary tree in `map[1..2n)`;
/// `map[0]` is unused. The element `a[i]` is stored as `value[i]` together with the maps on the
/// path from the leaf `n + i` up to the root, which have not been applied to it yet.
/// - `a[i] = act(map[1], ..., act(map[(n + i) / 2], act(map[n + i], value[i])) ...)`,
///   i.e. the maps are applied from the leaf to the root.
/// - The map at a node is newer than the maps at its descendants.
///
/// where `a` denotes the logical sequence of elements.
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
    /// Constructs a dual segment tree from a vector `v` of elements.
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

    /// Sets the element at index `i` to `act(f, a[i])`, where `a[i]` is the current element.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i` is out of bounds.
    pub fn apply(&mut self, i: usize, f: &M::Value) {
        assert!(
            i < self.len(),
            "index out of bounds: i={i}, len={}",
            self.len()
        );
        self.range_apply(i..=i, f);
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
    pub fn range_apply(&mut self, range: impl std::ops::RangeBounds<usize>, f: &M::Value) {
        let (mut l, mut r) = to_half_open(self.len(), range);
        assert!(
            l <= r,
            "left bound must be less than or equal to right bound: l={l}, r={r}"
        );
        assert!(r <= self.len(), "range out of bounds: range=[{l}, {r})");
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

    /// Returns the number of elements.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.value.len()
    }

    /// Returns `true` if the dual segment tree contains no elements.
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
            let k = i >> t;
            let id = self.monoid.id();
            let f = std::mem::replace(&mut self.map[k], id);
            self.map[k << 1] = self.monoid.op(&self.map[k << 1], &f);
            self.map[(k << 1) | 1] = self.monoid.op(&self.map[(k << 1) | 1], &f);
        }
    }
}

impl<T: Clone, M: Monoid, F: Action<T, M::Value>> DualSegmentTree<T, M, F> {
    /// Returns the element at index `i`, or `None` if `i` is out of bounds.
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
