/// A partial step function `f: K -> V` whose domain is a finite union of half-open intervals.
///
/// # Definition
/// A finite family of pairwise disjoint half-open intervals `[a_i, b_i)` of a totally ordered
/// set `K`, together with a value `v_i` in `V` for each interval. `f` is defined on the union
/// of the intervals and takes the constant value `v_i` on `[a_i, b_i)`. `m` is the number of
/// intervals.
///
/// Adjacent intervals with equal values are not merged automatically; see [`RangeMap::merge`].
///
/// # Invariants
/// - `self.0` sends `a_i` to `(b_i, v_i)` with `a_i < b_i`, and consecutive keys `a_i < a_j`
///   satisfy `b_i <= a_j`.
///
/// # Complexity
/// - Space: O(m)
pub struct RangeMap<K, V>(std::collections::BTreeMap<K, (K, V)>);

impl<K, V> RangeMap<K, V> {
    /// The map `f` defined nowhere.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn new() -> Self {
        Self(std::collections::BTreeMap::new())
    }

    /// The number `m` of intervals.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Whether `m = 0`, that is `f` is defined nowhere.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        self.0.len() == 0
    }
}

impl<K: Copy, V> RangeMap<K, V> {
    /// All the intervals in increasing order, as `(a, b, &v)`.
    ///
    /// # Complexity
    /// - Time: O(m) to exhaust
    /// - Space: O(1)
    pub fn iter(&self) -> impl DoubleEndedIterator<Item = (K, K, &V)> {
        self.0.iter().map(|(&a, (b, v))| (a, *b, v))
    }
}

impl<K: Ord + Copy, V: Clone> RangeMap<K, V> {
    /// The interval `[a, b)` containing `p` and its value `v`, as `(a, b, &v)`, or `None` if `f(p)`
    /// is undefined.
    ///
    /// # Complexity
    /// - Time: O(log m)
    /// - Space: O(1)
    pub fn get(&self, p: K) -> Option<(K, K, &V)> {
        let (&a, (b, v)) = self.0.range(..=p).next_back()?;
        (p < *b).then_some((a, *b, v))
    }

    /// Sets `f` to `v` on `[l, r)`, and returns the previous intervals clipped to `[l, r)`, in
    /// increasing order.
    ///
    /// # Complexity
    /// - Time: O((k + 1) log m), where `k` is the number of pieces returned, and amortized O(log m)
    /// - Space: O(k)
    ///
    /// # Panics
    /// Panics if `l >= r`.
    pub fn assign(&mut self, l: K, r: K, v: V) -> Vec<(K, K, V)> {
        let removed = self.remove(l, r);
        self.0.insert(l, (r, v));
        removed
    }

    /// Makes `f` undefined on `[l, r)`, and returns the previous intervals clipped to `[l, r)`, in
    /// increasing order.
    ///
    /// # Complexity
    /// - Time: O((k + 1) log m), where `k` is the number of intervals returned, and amortized
    ///   O(log m)
    /// - Space: O(k)
    ///
    /// # Panics
    /// Panics if `l >= r`.
    pub fn remove(&mut self, l: K, r: K) -> Vec<(K, K, V)> {
        assert!(l < r, "l must be less than r");
        self.split(l);
        self.split(r);
        let mut removed = vec![];
        while let Some((&a, _)) = self.0.range(l..r).next() {
            let (b, w) = self.0.remove(&a).unwrap();
            removed.push((a, b, w));
        }
        removed
    }

    /// The intervals intersecting `[l, r)`, unclipped, in increasing order, as `(a, b, &v)`.
    ///
    /// # Complexity
    /// - Time: O(log m + k) to exhaust, where `k` is the number of intervals yielded
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `l >= r`.
    pub fn range(&self, l: K, r: K) -> impl DoubleEndedIterator<Item = (K, K, &V)> {
        assert!(l < r, "l must be less than r");
        let head = self.0.range(..l).next_back().filter(|(_, (b, _))| l < *b);
        head.into_iter()
            .chain(self.0.range(l..r))
            .map(|(&a, (b, v))| (a, *b, v))
    }

    /// Makes `p` an endpoint of the intervals without changing `f`: if `a < p < b` for an interval
    /// `[a, b)`, replaces it by `[a, p)` and `[p, b)` with the same value, and otherwise does
    /// nothing.
    ///
    /// # Complexity
    /// - Time: O(log m)
    /// - Space: O(1)
    pub fn split(&mut self, p: K) {
        let Some((&a, entry)) = self.0.range_mut(..=p).next_back() else {
            return;
        };
        if !(a < p && p < entry.0) {
            return;
        }
        let b = std::mem::replace(&mut entry.0, p);
        let v = entry.1.clone();
        self.0.insert(p, (b, v));
    }
}

impl<K: Ord + Copy, V: Clone + PartialEq> RangeMap<K, V> {
    /// Merges the intervals forming the maximal interval containing `p` on which `f` is constant
    /// into one, and returns it, or `None` if `f(p)` is undefined.
    ///
    /// # Complexity
    /// - Time: O((k + 1) log m), where `k` is the number of intervals merged
    /// - Space: O(1)
    pub fn merge(&mut self, p: K) -> Option<(K, K)> {
        let (mut a, mut b, _) = self.get(p)?;
        while let Some((&c, &(d, ref w))) = self.0.range(..a).next_back() {
            if d != a || *w != self.0[&a].1 {
                break;
            }
            self.0.remove(&a);
            self.0.get_mut(&c).unwrap().0 = b;
            a = c;
        }
        while let Some((&c, &(d, ref w))) = self.0.range(b..).next() {
            if c != b || *w != self.0[&a].1 {
                break;
            }
            self.0.remove(&c);
            self.0.get_mut(&a).unwrap().0 = d;
            b = d;
        }
        Some((a, b))
    }
}

impl<K: Ord + Copy, V: Clone> Default for RangeMap<K, V> {
    fn default() -> Self {
        Self::new()
    }
}
