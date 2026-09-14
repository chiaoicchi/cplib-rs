use crate::algebra::{Bounded, Zero};
use crate::graph::max_flow::MaxFlow;

/// Minimization of a submodular quadratic function of binary variables, by a minimum cut.
///
/// # Definition
/// For variables `x_0, ..., x_{n-1}` in `{0, 1}`, the function is a sum of terms of two kinds:
/// a unary term `f_i(x_i)` given by its two values, and a binary term `f_{ij}(x_i, x_j)` given by
/// its four values, which must be submodular: `f(0, 0) + f(1, 1) <= f(0, 1) + f(1, 0)`. `solve`
/// returns the minimum of the sum and a minimizer.
///
/// Every such function is a constant plus a sum of non-negative terms of the forms `[x_i = 1] c`,
/// `[x_i = 0] c`, and `[x_i = 0, x_j = 1] c`, which are the capacities of the edges `s -> i`,
/// `i -> t`, and `i -> j` of a cut problem; the minimum is the constant plus the minimum cut, and
/// `x_i = 1` exactly when `i` lies on the sink side.
///
/// `T::max_value()` denotes `+∞` and forbids an assignment. A unary term may be infinite at one of
/// its values, and a binary term at `f(0, 1)` or `f(1, 0)`. The minimum is `+∞` if no assignment is
/// feasible. Finite costs must have a sum that does not overflow `T`.
///
/// To maximize, negate every cost.
///
/// # Invariants
/// - `base` is the accumulated constant, possibly `+∞`.
/// - `edges` holds `(u, v, c)` with `c > 0`, where `u`, `v` are variables or `s = n`, `t = n + 1`.
///
/// # Complexity
/// - Space: O(n + k), where `k` is the number of terms added
pub struct BinaryOptimization<T> {
    n: usize,
    base: T,
    edges: Vec<(usize, usize, T)>,
}

impl<T: Copy + Ord + Zero + Bounded + std::ops::Add<Output = T> + std::ops::Sub<Output = T>>
    BinaryOptimization<T>
{
    /// Constructs the zero function of `n` variables.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn new(n: usize) -> Self {
        Self {
            n,
            base: T::zero(),
            edges: Vec::new(),
        }
    }

    /// Adds the unary term with value `x0` at `x_i = 0` and `x1` at `x_i = 1`.
    ///
    /// # Complexity
    /// - Time: O(1) amortized
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i >= n` or both values are `+∞`.
    pub fn add_unary(&mut self, i: usize, x0: T, x1: T) {
        assert!(i < self.n, "index out of bounds: i={i}, n={}", self.n);
        assert!(
            x0 != T::max_value() || x1 != T::max_value(),
            "a unary term must be finite at some value"
        );
        if x0 <= x1 {
            self.base = add(self.base, x0);
            self.add_edge(self.n, i, sub(x1, x0));
        } else {
            self.base = add(self.base, x1);
            self.add_edge(i, self.n + 1, sub(x0, x1));
        }
    }

    /// Adds the binary term with value `x_ab` at `(x_i, x_j) = (a, b)`.
    ///
    /// # Complexity
    /// - Time: O(1) amortized
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i >= n`, `j >= n`, `i == j`, `x00` or `x11` is `+∞`, or the term
    /// is not submodular.
    pub fn add_binary(&mut self, i: usize, j: usize, x00: T, x01: T, x10: T, x11: T) {
        assert!(i < self.n, "index out of bounds: i={i}, n={}", self.n);
        assert!(j < self.n, "index out of bounds: j={j}, n={}", self.n);
        assert!(i != j, "variables must differ: i={i}");
        assert!(
            x00 != T::max_value() && x11 != T::max_value(),
            "a binary term must be finite at (0, 0) and (1, 1)"
        );
        assert!(add(x00, x11) <= add(x01, x10), "the term is not submodular");
        if x01 == T::max_value() && x10 == T::max_value() {
            self.add_unary(i, x00, x11);
            self.add_edge(i, j, T::max_value());
            self.add_edge(j, i, T::max_value());
        } else if x10 == T::max_value() {
            self.add_binary(j, i, x00, x10, x01, x11);
        } else {
            self.add_unary(i, x00, x10);
            self.add_unary(j, T::zero(), sub(x11, x10));
            self.add_edge(i, j, sub(add(x01, x10), add(x00, x11)));
        }
    }

    /// Returns the minimum value and a minimizer.
    ///
    /// # Complexity
    /// - Time: O(n^2 k)
    /// - Space: O(n + k)
    pub fn solve(self) -> (T, Vec<bool>) {
        let mut max_flow = MaxFlow::from_edges(self.n + 2, &self.edges);
        let value = add(self.base, max_flow.flow(self.n, self.n + 1));
        (
            value,
            max_flow
                .reachable(self.n)
                .iter()
                .take(self.n)
                .map(|x| !x)
                .collect(),
        )
    }

    fn add_edge(&mut self, u: usize, v: usize, c: T) {
        if c > T::zero() {
            self.edges.push((u, v, c));
        }
    }
}

/// `a + b`, where `T::max_value()` is `+∞`.
fn add<T: Copy + Ord + Bounded + std::ops::Add<Output = T>>(a: T, b: T) -> T {
    if a == T::max_value() || b == T::max_value() {
        T::max_value()
    } else {
        a + b
    }
}

/// `a - b` for `b` finite, where `T::max_value()` is `+∞`.
fn sub<T: Copy + Ord + Bounded + std::ops::Sub<Output = T>>(a: T, b: T) -> T {
    assert!(b != T::max_value());
    if a == T::max_value() {
        T::max_value()
    } else {
        a - b
    }
}
