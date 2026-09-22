use crate::algebra::Group;

/// A partition of `[0, n)` into disjoint sets with potentials in a group, with unions of its sets
/// by constraints on the potentials.
///
/// # Definition
/// `S(x)` is the set containing `x`. A potential on a set `S` is a map `φ: S -> G`, and the
/// potential to `y` from `x` is `op(inv(φ(x)), φ(y))`. Each set carries a potential determined, up
/// to [`unite`](Self::unite).
///
/// # Invariants
/// - `value` is a rooted forest on `[0, n)` whose trees are the sets: `value[i]` is the parent of
///   `i` if `value[i] >= 0`, and otherwise `i` is a root and `-value[i]` is the size of its tree.
/// - `potential[i]` is the potential to `i` from its parent if `i` is not a root, and unused
///   otherwise.
/// - `count` is the number of sets.
///
/// # Complexity
/// - Space: O(n)
pub struct PotentialDsu<G: Group> {
    group: G,
    value: Box<[i32]>,
    potential: Box<[G::Value]>,
    count: usize,
}

impl<G: Group> PotentialDsu<G> {
    /// The partition of `[0, n)` into singletons.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `n >= 2^31`.
    pub fn new(group: G, n: usize) -> Self {
        assert!(n < 1 << 31, "n must be less than 2^31: n={n}");
        Self {
            potential: (0..n).map(|_| group.id()).collect(),
            group,
            value: vec![-1; n].into(),
            count: n,
        }
    }

    /// The representative `r` of `S(x)`, as in [`Dsu::root`](crate::collections::dsu::Dsu::root),
    /// and the potential to `x` from `r`.
    ///
    /// # Complexity
    /// - Time: amortized O(α(n)), where `α` is the inverse Ackermann function
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x >= n`.
    pub fn root(&mut self, mut x: usize) -> (usize, G::Value) {
        assert!(
            x < self.len(),
            "index out of bounds: x={x}, len={}",
            self.len()
        );
        let mut p = self.group.id();
        while self.value[x] >= 0 {
            let px = self.value[x] as usize;
            p = self.group.op(&self.potential[x], &p);
            if self.value[px] >= 0 {
                self.value[x] = self.value[px];
                self.potential[x] = self.group.op(&self.potential[px], &self.potential[x]);
            }
            x = px;
        }
        (x, p)
    }

    /// Imposes that the potential to `y` from `x` is `p`, and returns whether this is consistent
    /// with the constraints so far.
    ///
    /// # Definition
    /// If `S(x) != S(y)`, unites the two sets and returns `true`. Otherwise changes nothing and
    /// returns whether the potential to `y` from `x` is `p`.
    ///
    /// # Complexity
    /// - Time: amortized O(α(n))
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x >= n` or `y >= n`.
    pub fn unite(&mut self, x: usize, y: usize, p: &G::Value) -> bool
    where
        G::Value: PartialEq,
    {
        let (mut rx, potx) = self.root(x);
        let (mut ry, poty) = self.root(y);
        if rx == ry {
            return self.group.op(&potx, p) == poty;
        }
        let mut p = self
            .group
            .op(&self.group.op(&potx, p), &self.group.inv(&poty));
        if self.value[rx] > self.value[ry] {
            std::mem::swap(&mut rx, &mut ry);
            p = self.group.inv(&p);
        }
        self.value[rx] += self.value[ry];
        self.value[ry] = rx as i32;
        self.potential[ry] = p;
        self.count -= 1;
        true
    }

    /// The potential to `y` from `x`, or `None` if `S(x) != S(y)`.
    ///
    /// # Complexity
    /// - Time: amortized O(α(n))
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x >= n` or `y >= n`.
    pub fn potential(&mut self, x: usize, y: usize) -> Option<G::Value> {
        let (rx, potx) = self.root(x);
        let (ry, poty) = self.root(y);
        if rx == ry {
            Some(self.group.op(&self.group.inv(&potx), &poty))
        } else {
            None
        }
    }

    /// The size `|S(x)|`.
    ///
    /// # Complexity
    /// - Time: amortized O(α(n))
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x >= n`.
    pub fn set_size(&mut self, x: usize) -> usize {
        -self.value[self.root(x).0] as usize
    }

    /// The number of sets.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn num_sets(&self) -> usize {
        self.count
    }

    /// The number `n` of elements.
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
}
