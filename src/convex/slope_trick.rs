/// A convex piecewise linear function `f` with integer slopes (slope trick).
///
/// # Definition
/// `f: Z -> Z` is convex and piecewise linear with integer slopes, and is bounded below.
///
/// # Invariants
/// - `n` is the number of breakpoints with multiplicity. `left` holds `l - shift_left` for the
///   breakpoints `l` at which the slope decreases goint left, and `right` holds `r - shift_right`
///   for the breakpoints `r` at which it increases going right, so that
///   `f(x) = min + Σ_l max(0, l - x) + Σ_r min(0, x - r)`.
/// - `min` is the minimum of `f`.
///
/// # Complexity
/// - Space: O(n)
pub struct SlopeTrick {
    left: std::collections::BinaryHeap<i64>,
    right: std::collections::BinaryHeap<std::cmp::Reverse<i64>>,
    shift_left: i64,
    shift_right: i64,
    min: i64,
}
impl Default for SlopeTrick {
    fn default() -> Self {
        Self::new()
    }
}

impl SlopeTrick {
    /// The zero function.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn new() -> Self {
        Self {
            left: std::collections::BinaryHeap::new(),
            right: std::collections::BinaryHeap::new(),
            shift_left: 0,
            shift_right: 0,
            min: 0,
        }
    }

    /// The minimum of `f`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn min(&self) -> i64 {
        self.min
    }

    /// The interval of which `f` attains its minimum, with `None` for an unbounded end.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn argmin(&self) -> (Option<i64>, Option<i64>) {
        (self.top_left(), self.top_right())
    }

    /// Adds the constant `c` to `f`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn add_const(&mut self, c: i64) {
        self.min += c;
    }

    /// Adds `max(0, a - x)` to `f`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    pub fn add_left(&mut self, a: i64) {
        if let Some(r) = self.top_right() {
            self.min += (a - r).max(0);
        }
        self.push_right(a);
        let x = self.pop_right().unwrap();
        self.push_left(x);
    }

    /// Adds `max(0, x - a)` to `f`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    pub fn add_right(&mut self, a: i64) {
        if let Some(l) = self.top_left() {
            self.min += (l - a).max(0);
        }
        self.push_left(a);
        let x = self.pop_left().unwrap();
        self.push_right(x);
    }

    /// Adds `|x - a|` to `f`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    pub fn add_abs(&mut self, a: i64) {
        self.add_left(a);
        self.add_right(a);
    }

    /// Adds `ax + b` to `f`.
    ///
    /// # Complexity
    /// - Time: O(|a| log n)
    /// - Space: O(|a|)
    ///
    /// # Panics
    /// Panics if `a > num_left()` or `-a > num_right()`, that is if `f + ax + b` is unbounded
    /// below.
    pub fn add_linear(&mut self, a: i64, b: i64) {
        assert!(
            a <= self.num_left() as i64 && -a <= self.num_right() as i64,
            "f + ax + b must be bounded below: a={a}, num_left={}, num_right={}",
            self.num_left(),
            self.num_right(),
        );
        self.min += b;
        for _ in 0..a.max(0) {
            let x = self.pop_left().unwrap();
            self.min += x;
            self.push_right(x);
        }
        for _ in 0..(-a).max(0) {
            let x = self.pop_right().unwrap();
            self.min -= x;
            self.push_left(x);
        }
    }

    /// Replaces `f` by `g(x) = min_{y<=x} f(y)`, which is non-increasing.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(1)
    pub fn flatten_right(&mut self) {
        self.right.clear();
    }

    /// Replaces `f` by `g(x) = max_{y>=x} f(y)`, which is non-decreasing.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(1)
    pub fn flatten_left(&mut self) {
        self.left.clear();
    }

    /// Replaces `f` by `g(x) = f(x - a)`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn shift(&mut self, a: i64) {
        self.shift_left += a;
        self.shift_right += a;
    }

    /// Replaces `f` by `g(x) = min_{x+lo<=y<=x+hi} f(y)`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `lo > hi`.
    pub fn window_min(&mut self, lo: i64, hi: i64) {
        assert!(lo <= hi, "lo must be at most hi: lo={lo}, hi={hi}");
        self.shift_left -= hi;
        self.shift_right -= lo;
    }

    /// The value `f(x)`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(1)
    pub fn eval(&self, x: i64) -> i64 {
        let mut y = self.min;
        for &l in &self.left {
            y += (l + self.shift_left - x).max(0);
        }
        for &std::cmp::Reverse(r) in &self.right {
            y += (x - (r + self.shift_right)).max(0);
        }
        y
    }

    /// The number of breakpoints at which the slope decreases going left, with multiplicity.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn num_left(&self) -> usize {
        self.left.len()
    }

    /// The number of breakpoints at which the slope increases going right, with multiplicity.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn num_right(&self) -> usize {
        self.right.len()
    }

    /// The number `n` of breakpoints, with multiplicity.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn num_breakpoints(&self) -> usize {
        self.num_left() + self.num_right()
    }

    /// The greatest breakpoint at which the slope decreses going left, or `None` if there is none.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    fn top_left(&self) -> Option<i64> {
        self.left.peek().map(|&x| x + self.shift_left)
    }

    /// The least breakpoint at which the slope increases going right, or `None` if there is none.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    fn top_right(&self) -> Option<i64> {
        self.right
            .peek()
            .map(|&std::cmp::Reverse(x)| x + self.shift_right)
    }

    /// Adds `x` to the breakpoints at which the slope decreases goint left.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: amortized O(1)
    fn push_left(&mut self, x: i64) {
        self.left.push(x - self.shift_left)
    }

    /// Adds `x` to the breakpoints at which the slope increases goint right.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: amortized O(1)
    fn push_right(&mut self, x: i64) {
        self.right.push(std::cmp::Reverse(x - self.shift_right));
    }

    /// Removes the greatest breakpoint at which the slope decreses going left and returns it, or
    /// `None` if there is none.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    fn pop_left(&mut self) -> Option<i64> {
        self.left.pop().map(|x| x + self.shift_left)
    }

    /// Removes the least breakpoint at which the slope increases going right and returns it, or
    /// `None` if there is none.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    fn pop_right(&mut self) -> Option<i64> {
        self.right
            .pop()
            .map(|std::cmp::Reverse(x)| x + self.shift_right)
    }
}
