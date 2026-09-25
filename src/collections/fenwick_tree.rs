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
    pub fn max_right(&self, l: usize, mut pred: impl FnMut(&G::Value) -> bool) -> usize {
        let n = self.len();
        assert!(l <= n, "index out of bounds: l={l}, n={n}");
        assert!(pred(&self.monoid.id()), "pred(id()) must be true");
        let inv = self.monoid.inv(&self.prefix_fold(l));
        let mut r = 0;
        let mut acc = self.monoid.id();
        let mut k = (n + 1).next_power_of_two() >> 1;
        while 0 < k {
            if r + k <= n {
                let next = self.monoid.op(&acc, &self.value[r + k]);
                if r + k <= l || pred(&self.monoid.op(&inv, &next)) {
                    acc = next;
                    r += k;
                }
            }
            k >>= 1;
        }
        r
    }

    /// An `l` in `[0, r]` at which `pred(fold[l, r))` turns from `true` to `false` as `l` decreases.
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
    pub fn min_left(&self, r: usize, mut pred: impl FnMut(&G::Value) -> bool) -> usize {
        let n = self.len();
        assert!(r <= n, "index out of bounds: r={r}, len={n}");
        assert!(pred(&self.monoid.id()), "pred(id()) must be true");
        let whole = self.prefix_fold(r);
        if pred(&whole) {
            return 0;
        }
        let mut l = 0;
        let mut acc = self.monoid.id();
        let mut k = (r + 1).next_power_of_two() >> 1;
        while 0 < k {
            if l + k < r {
                let next = self.monoid.op(&acc, &self.value[l + k]);
                if !pred(&self.monoid.op(&self.monoid.inv(&next), &whole)) {
                    acc = next;
                    l += k;
                }
            }
            k >>= 1;
        }
        l + 1
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
