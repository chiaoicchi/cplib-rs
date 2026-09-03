/// A disjoint set union data structure.
///
/// # Invariants
/// The sets form a partition of `{0, 1, ..., n - 1}`, represented by a rooted forest
/// on `{0, 1, ..., n - 1}` stored in the internal array `value`.
/// - If `value[i] < 0`, then `i` is a root and `-value[i]` is the size of its set.
/// - If `value[i] >= 0`, then `value[i]` is the parent of `i`.
///
/// Two elements belong to the same set if and only if they belong to the same tree, and the root of
/// the tree is the representative of the set.
///
/// # Complexity
/// - Space: O(n)
pub struct Dsu {
    value: Box<[i32]>,
    count: usize,
}

impl Dsu {
    /// Constructs a disjoint set union with `n` elements, each in its own set.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `n` is greater than or equal to `2^31`.
    pub fn new(n: usize) -> Self {
        assert!(n < 1 << 31, "n must be less than 2^31: n={n}");
        Self {
            value: vec![-1; n].into_boxed_slice(),
            count: n,
        }
    }

    /// Returns the representative of the set containing `x`.
    ///
    /// # Complexity
    /// - Time: amortized O(α(n)), where `α` is the inverse Ackermann function
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x` is out of bounds.
    pub fn root(&mut self, mut x: usize) -> usize {
        assert!(
            x < self.len(),
            "index out of bounds: x={x}, len={}",
            self.len()
        );
        while self.value[x] >= 0 {
            let px = self.value[x] as usize;
            if self.value[px] >= 0 {
                self.value[x] = self.value[px];
            }
            x = px;
        }
        x
    }

    /// Unites the sets containing `x` and `y`, and returns `true` if they were different sets.
    ///
    /// # Complexity
    /// - Time: amortized O(α(n))
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x` is out of bounds.
    /// Panics if `y` is out of bounds.
    pub fn unite(&mut self, x: usize, y: usize) -> bool {
        let mut rx = self.root(x);
        let mut ry = self.root(y);
        if rx == ry {
            return false;
        }
        if self.value[rx] > self.value[ry] {
            std::mem::swap(&mut rx, &mut ry);
        }
        self.value[rx] += self.value[ry];
        self.value[ry] = rx as i32;
        self.count -= 1;
        true
    }

    /// Returns `true` if `x` and `y` belong to the same set.
    ///
    /// # Complexity
    /// - Time: amortized O(α(n))
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x` is out of bounds.
    /// Panics if `y` is out of bounds.
    pub fn is_same(&mut self, x: usize, y: usize) -> bool {
        self.root(x) == self.root(y)
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
        -self.value[self.root(x)] as usize
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

    /// Returns `true` if the disjoint set union contains no elements.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}
