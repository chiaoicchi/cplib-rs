/// A partition of `[0, n)` into disjoint sets, with unions of its sets.
///
/// # Definition
/// `S(x)` is the set containing `x`.
///
/// # Invariants
/// - `value` is a rooted forest on `[0, n)` whose trees are the sets: `value[i]` is the parent of
///   `i` if `value[i] >= 0`, and otherwise `i` is a root and `-value[i]` is the size of its tree.
/// - `count` is the number of sets.
///
/// # Complexity
/// - Space: O(n)
pub struct Dsu {
    value: Box<[i32]>,
    count: usize,
}

impl Dsu {
    /// The partition of `[0, n)` into singletons.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `n >= 2^31`.
    pub fn new(n: usize) -> Self {
        assert!(n < 1 << 31, "n must be less than 2^31: n={n}");
        Self {
            value: vec![-1; n].into_boxed_slice(),
            count: n,
        }
    }

    /// The representative of `S(x)`, an element of `S(x)` shared by all its elements, which stays
    /// the same until `S(x)` is united with another set.
    ///
    /// # Complexity
    /// - Time: amortized O(α(n)), where `α` is the inverse Ackermann function
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x >= n`.
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

    /// Unites `S(x)` and `S(y)`, and returns whether they were different.
    ///
    /// # Complexity
    /// - Time: amortized O(α(n))
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x >= n` or `y >= n`.
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

    /// Whether `S(x) = S(y)`.
    ///
    /// # Complexity
    /// - Time: amortized O(α(n))
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x >= n` or `y >= n`.
    pub fn is_same(&mut self, x: usize, y: usize) -> bool {
        self.root(x) == self.root(y)
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
        -self.value[self.root(x)] as usize
    }

    /// The numbers of sets.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn num_sets(&self) -> usize {
        self.count
    }

    /// Whether `n = 0`.
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
