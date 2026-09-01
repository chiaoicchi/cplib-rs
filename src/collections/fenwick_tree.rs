use crate::algebra::{Commutative, Group};
use crate::range::to_half_open;

/// A Fenwick tree data structure.
///
/// # Invariants
/// The elements are stored in the internal array `value` as a 1-indexed array `value[1..=n]`;
/// `value[0]` is unused. Let `a` be a 1-indexed array `a[1..=n]`.
/// - `value[i] = op(a[i - lsb(i) + 1], a[i - lsb(i) + 2], ..., a[i])`,
///   where `lsb` (least significant bit) is the largest power of two that divides `i`.
///
/// # Complexity
/// - Space: O(n)
pub struct FenwickTree<G: Group + Commutative> {
    group: G,
    value: Vec<G::Value>,
}

impl<G: Group + Commutative> FenwickTree<G> {
    /// Constructs a Fenwick tree with `n` elements, all initialized to `id()`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    pub fn new(group: G, n: usize) -> Self {
        Self {
            value: (0..=n).map(|_| group.id()).collect(),
            group,
        }
    }

    /// Constructs a Fenwick tree from a vector `v` of elements.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    pub fn from_vec(group: G, mut v: Vec<G::Value>) -> Self {
        let n = v.len();
        v.insert(0, group.id());
        for i in 1..n {
            let lsb = lsb(i);
            if i + lsb <= n {
                v[i + lsb] = group.op(&v[i], &v[i + lsb]);
            }
        }
        Self { group, value: v }
    }

    /// Sets the element at index `i` to `x`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i` is out of bounds.
    pub fn set(&mut self, i: usize, x: &G::Value) {
        assert!(
            i < self.len(),
            "index out of bounds: i={i}, len={}",
            self.len()
        );
        let a = self.get(i).unwrap();
        let x = self.group.op(&self.group.inv(&a), x);
        self.op_assign(i, &x);
    }

    /// Sets the element at index `i` to `op(a[i], x)`, where `a[i]` is the current element.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i` is out of bounds.
    pub fn op_assign(&mut self, mut i: usize, x: &G::Value) {
        assert!(
            i < self.len(),
            "index out of bounds: i={i}, len={}",
            self.len()
        );
        i += 1;
        while i < self.value.len() {
            self.value[i] = self.group.op(&self.value[i], x);
            i += lsb(i);
        }
    }

    /// Appends the element `x` to the back of the Fenwick tree.
    ///
    /// # Complexity
    /// - Time: worst O(log n), average: O(1)
    /// - Space: worst O(n), average: O(1)
    pub fn push(&mut self, mut x: G::Value) {
        let n = self.value.len();
        let mut k = 1;
        let lsb = lsb(n);
        while k < lsb {
            x = self.group.op(&self.value[n - k], &x);
            k <<= 1;
        }
        self.value.push(x);
    }

    /// Returns the value at index `i`, or `None` if `i` is out of bounds.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    pub fn get(&self, i: usize) -> Option<G::Value> {
        if self.len() <= i {
            None
        } else {
            Some(self.fold(i..=i))
        }
    }

    /// Folds the elements in `0..r`.
    ///
    /// # Panics
    /// Panics if `r` is out of bounds.
    fn prefix_fold(&self, mut r: usize) -> G::Value {
        assert!(
            r <= self.len(),
            "index out of bounds: r={r}, len={}",
            self.len()
        );
        let mut x = self.group.id();
        while 0 < r {
            x = self.group.op(&self.value[r], &x);
            r -= lsb(r);
        }
        x
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
    pub fn fold(&self, range: impl std::ops::RangeBounds<usize>) -> G::Value {
        let (l, r) = to_half_open(self.len(), range);
        assert!(
            l <= r,
            "left bound must be less than or equal to right bound: l={}, r={}",
            l,
            r
        );
        assert!(r <= self.len(), "range out of bounds: range=[{l}, {r})");
        self.group
            .op(&self.group.inv(&self.prefix_fold(l)), &self.prefix_fold(r))
    }

    /// Returns the number of elements.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.value.len() - 1
    }

    /// Returns `true` if the Fenwick tree contains no elements.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Returns the least significant bit of `i`, i.e. the largest power of two that divides `i`.
#[inline]
fn lsb(i: usize) -> usize {
    i & i.wrapping_neg()
}
