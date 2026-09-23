use crate::algebra::Semigroup;

/// A deque of a semigroup, with the fold of its elements (sliding window aggregation).
///
/// # Definition
/// `n` is the number of elements. For the elements `x_0, ..., x_{n-1}` from the front, the fold is
/// `x_0 x_1 ... x_{n-1}`, the product under `op` in this order.
///
/// # Invariants
/// - `front_raw` holds the elements at the front in reverse order, and `back_raw` those at the back
///   in order, so that the queue is `front_raw` reversed followed by `back_raw`.
/// - `front_fold[i]` is the fold of `front_raw[..=i]` reversed, and `back_fold[i]` is the fold of
///   `back_raw[..=i]`.
///
/// # Complexity
/// - Space: O(n)
pub struct FoldableDeque<S: Semigroup> {
    semigroup: S,
    front_raw: Vec<S::Value>,
    front_fold: Vec<S::Value>,
    back_raw: Vec<S::Value>,
    back_fold: Vec<S::Value>,
}

impl<S: Semigroup<Value: Clone>> FoldableDeque<S> {
    /// The empty deque.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn new(semigroup: S) -> Self {
        Self {
            semigroup,
            front_raw: Vec::new(),
            front_fold: Vec::new(),
            back_raw: Vec::new(),
            back_fold: Vec::new(),
        }
    }

    /// Prepends `x` to the front.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: amortized O(1)
    pub fn push_front(&mut self, x: S::Value) {
        let fold = match self.front_fold.last() {
            Some(y) => self.semigroup.op(&x, y),
            None => x.clone(),
        };
        self.front_raw.push(x);
        self.front_fold.push(fold);
    }

    /// Appends `x` to the back.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: amortized O(1)
    pub fn push_back(&mut self, x: S::Value) {
        let fold = match self.back_fold.last() {
            Some(y) => self.semigroup.op(y, &x),
            None => x.clone(),
        };
        self.back_raw.push(x);
        self.back_fold.push(fold);
    }

    /// Removes the element at the front and returns it, or `None` if the deque is empty.
    ///
    /// # Complexity
    /// - Time: amortized O(1)
    /// - Space: O(1)
    pub fn pop_front(&mut self) -> Option<S::Value> {
        if self.front_raw.is_empty() {
            let half = self.back_raw.len().div_ceil(2);
            let rest = self.back_raw.split_off(half);
            self.front_raw = self.back_raw.split_off(0);
            self.front_raw.reverse();
            self.back_raw = rest;
            self.rebuild();
        }
        self.front_fold.pop();
        self.front_raw.pop()
    }

    /// Removes the element at the back and returns it, or `None` if the deque is empty.
    ///
    /// # Complexity
    /// - Time: amortized O(1)
    /// - Space: O(1)
    pub fn pop_back(&mut self) -> Option<S::Value> {
        if self.back_raw.is_empty() {
            let half = self.front_raw.len().div_ceil(2);
            let rest = self.front_raw.split_off(half);
            self.back_raw = self.front_raw.split_off(0);
            self.back_raw.reverse();
            self.front_raw = rest;
            self.rebuild();
        }
        self.back_fold.pop();
        self.back_raw.pop()
    }

    /// Recomputes `front_fold` and `back_fold` from `front_raw` and `back_raw`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(1)
    fn rebuild(&mut self) {
        self.front_fold.clear();
        for x in &self.front_raw {
            let fold = match self.front_fold.last() {
                Some(y) => self.semigroup.op(x, y),
                None => x.clone(),
            };
            self.front_fold.push(fold);
        }
        self.back_fold.clear();
        for x in &self.back_raw {
            let fold = match self.back_fold.last() {
                Some(y) => self.semigroup.op(y, x),
                None => x.clone(),
            };
            self.back_fold.push(fold);
        }
    }

    /// The element at the front, or `None` if the deque is empty.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn front(&self) -> Option<&S::Value> {
        self.front_raw.last().or(self.back_raw.first())
    }

    /// The element at the back, or `None` if the deque is empty.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn back(&self) -> Option<&S::Value> {
        self.back_raw.last().or(self.front_raw.first())
    }

    /// The fold of the elements, or `None` if the deque is empty.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn fold(&self) -> Option<S::Value> {
        match (self.front_fold.last(), self.back_fold.last()) {
            (Some(x), Some(y)) => Some(self.semigroup.op(x, y)),
            (Some(x), None) => Some(x.clone()),
            (None, y) => y.cloned(),
        }
    }

    /// The number `n` of elements.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.front_raw.len() + self.back_raw.len()
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
