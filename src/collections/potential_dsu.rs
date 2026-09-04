use crate::algebra::Group;

/// A potential disjoint set union data structure.
///
/// # Definition
/// A potential on a set `S` is a map `φ: S -> G`.
/// The potential to `y` from `x` is `op(inv(φ(x)), φ(y))`; it is unchanged when `φ` is replaced by
/// `g.φ` for any `g` in `G`, so a potential is determined only up to left translation.
/// This structure maintains, on each set of the partition, a potential determined up to left
/// translation by the constraints given to [`unite`](Self::unite).
///
/// # Invariants
/// The sets form a partition of `{0, 1, ..., n - 1}`, represented by a rooted forest
/// on `{0, 1, ..., n - 1}` stored in the internal array `value`, as in
/// [`Dsu`](crate::collections::dsu::Dsu).
/// - If `value[i] < 0`, then `i` is a root and `-value[i]` is the size of its set.
/// - If `value[i] >= 0`, then `value[i]` is the parent of `i`, and `potential[i]` is the potential
///   to `i` from `value[i]`.
///
/// Hence the potential to `i` from the root of its tree is the product of `potential` along the
/// path from the root to `i`, with the root side on the left.
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
    /// Constructs a potential disjoint set union with `n` elements, each in its own set.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `n` is greater than or equal to `2^31`.
    pub fn new(group: G, n: usize) -> Self {
        assert!(n < 1 << 31, "n must be less than 2^31: n={n}");
        Self {
            potential: (0..n).map(|_| group.id()).collect(),
            group,
            value: vec![-1; n].into(),
            count: n,
        }
    }

    /// Returns the representative `r` of the set containing `x` and the potential to `x` from `r`.
    ///
    /// # Complexity
    /// - Time: amortized O(α(n)), where `α` is the inverse Ackermann function
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x` is out of bounds.
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

    /// Imposes the constraint that the potential to `y` from `x` is `p`.
    ///
    /// If `x` and `y` are in different sets, unites them and returns `true`.
    /// Otherwise leaves the partition unchanged and returns whether the constraint is consistent
    /// with the existing ones, i.e. whether `potential(x, y) == Some(p)`.
    ///
    /// # Complexity
    /// - Time: amortized O(α(n))
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x` or `y` is out of bounds.
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

    /// Returns the potential to `y` from `x`, or `None` if they are in different sets.
    ///
    /// # Complexity
    /// - Time: amortized O(α(n))
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x` or `y` is out of bounds.
    pub fn potential(&mut self, x: usize, y: usize) -> Option<G::Value> {
        let (rx, potx) = self.root(x);
        let (ry, poty) = self.root(y);
        if rx == ry {
            Some(self.group.op(&self.group.inv(&potx), &poty))
        } else {
            None
        }
    }

    /// Returns the size of the set containing `x`.
    ///
    /// # Complexity
    /// - Time: amortized O(α(n))
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x` is out of bounds.
    pub fn set_size(&mut self, x: usize) -> usize {
        -self.value[self.root(x).0] as usize
    }

    /// Returns the number of disjoint sets.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn num_sets(&self) -> usize {
        self.count
    }

    /// Returns the number of elements.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.value.len()
    }

    /// Returns `true` if the potential disjoint set union contains no elements.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}
