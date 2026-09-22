use crate::algebra::{Commutative, Group, Monoid};
use crate::range::to_half_open;

/// A sequence `a` of a monoid, with the folds of its prefixes.
///
/// # Definition
/// `n` is the length of `a`. `fold[l, r) = a[l] a[l + 1] ... a[r - 1]`, the product under `op` in
/// this order, which is `id()` for `l = r`.
///
/// # Invariants
/// - `value[i] = fold[i - lsb(i), i)` for `i` in `[1, n]`. `value[0]` is unused.
///
/// # Complexity
/// - Space: O(n)
pub struct FenwickTree<M: Monoid> {
    monoid: M,
    value: Vec<M::Value>,
}

impl<M: Monoid> FenwickTree<M> {
    /// The sequence of `n` copies of `id()`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    pub fn new(monoid: M, n: usize) -> Self {
        Self {
            value: (0..=n).map(|_| monoid.id()).collect(),
            monoid,
        }
    }

    /// The Fenwick tree of the sequence `v`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    pub fn from_vec(monoid: M, mut v: Vec<M::Value>) -> Self {
        v.insert(0, monoid.id());
        for j in 1..v.len() {
            let mut k = 1;
            while k < lsb(j) {
                v[j] = monoid.op(&v[j - k], &v[j]);
                k <<= 1;
            }
        }
        Self { monoid, value: v }
    }

    /// Appends `x` to `a`.
    ///
    /// # Complexity
    /// - Time: O(log n), amortized O(1)
    /// - Space: amortized O(1)
    pub fn push(&mut self, mut x: M::Value) {
        let n = self.value.len();
        let mut k = 1;
        let lsb = lsb(n);
        while k < lsb {
            x = self.monoid.op(&self.value[n - k], &x);
            k <<= 1;
        }
        self.value.push(x);
    }

    /// The fold `fold[0, r)`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `r > n`.
    pub fn prefix_fold(&self, mut r: usize) -> M::Value {
        assert!(
            r <= self.len(),
            "index out of bounds: r={r}, len={}",
            self.len()
        );
        let mut x = self.monoid.id();
        while 0 < r {
            x = self.monoid.op(&self.value[r], &x);
            r -= lsb(r);
        }
        x
    }

    /// The length `n`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.value.len() - 1
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

impl<M: Monoid + Commutative> FenwickTree<M> {
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
        i += 1;
        while i < self.value.len() {
            self.value[i] = self.monoid.op(&self.value[i], x);
            i += lsb(i);
        }
    }
}

impl<G: Group> FenwickTree<G> {
    /// The element `a[i]`, or `None` if `i >= n`.
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

    /// The fold `fold[l, r)` of `range = [l, r)`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `l > r` or `r > n`.
    pub fn fold(&self, range: impl std::ops::RangeBounds<usize>) -> G::Value {
        let (l, r) = to_half_open(self.len(), range);
        self.monoid
            .op(&self.monoid.inv(&self.prefix_fold(l)), &self.prefix_fold(r))
    }
}

impl<G: Group + Commutative> FenwickTree<G> {
    /// Sets `a[i]` to `x`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i >= n`.
    pub fn set(&mut self, i: usize, x: &G::Value) {
        assert!(
            i < self.len(),
            "index out of bounds: i={i}, len={}",
            self.len()
        );
        let a = self.get(i).unwrap();
        let x = self.monoid.op(&self.monoid.inv(&a), x);
        self.op_assign(i, &x);
    }
}

/// The largest power of two dividing `i`, and `0` for `i = 0`.
///
/// # Complexity
/// - Time: O(1)
/// - Space: O(1)
#[inline]
fn lsb(i: usize) -> usize {
    i & i.wrapping_neg()
}
