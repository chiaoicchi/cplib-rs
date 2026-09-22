use crate::algebra::{Bounded, Zero};
use crate::graph::max_flow::MaxFlow;

/// Minimization of a submodular quadratic function of binary variables, by a minimum cut.
///
/// # Definition
/// The function of `x_0, ..., x_{n-1}` in `{0, 1}` is a sum of unary terms `f_i(x_i)`, given by
/// their two values, and binary terms `f_{ij}(x_i, x_j)`, given by their four values, which are
/// submodular: `f(0, 0) + f(1, 1) <= f(0, 1) + f(1, 0)`. `k` is the number of terms.
///
/// `T::max_value()` is `+∞`, and forbids the assignments at which it is taken. A unary term may be
/// `+∞` at opne of its values, and a binary term at `f(0, 1)` or `f(0, 1)`.
///
/// # Contract
/// `4 Σ |v| < T::max_value()`, the sum over the finite values `v` of the terms.
///
/// # Invariants
/// - The function is `base` plus the value of the cut: the sum of `c` over the `(u, v, c)` in
///   `edges` with `u` on the source side and `v` on the sink side, where `s = n` is on the source
///   side, `t = n + 1` on the sink side, and `i` on the sink side iff `x_i = 1`.
/// - `c > 0` for every `(u, v, c)` in `edges`.
///
/// # Complexity
/// - Space: O(n + k)
pub struct BinaryOptimization<T> {
    n: usize,
    base: T,
    edges: Vec<(usize, usize, T)>,
}

impl<T: Copy + Ord + Zero + Bounded + std::ops::Add<Output = T> + std::ops::Sub<Output = T>>
    BinaryOptimization<T>
{
    /// The zero function of `n` variables.
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

    /// Adds the unary term of `x_i` with the value `x0` at `0` and `x1` and `1`.
    ///
    /// # Complexity
    /// - Time: amortized O(1)
    /// - Space: amortized O(1)
    ///
    /// # Panics
    /// Panics if `i >= n` or `x0 = x1 = T::max_value()`.
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

    /// Adds the binary term of `(x_i, x_j)` with the value `x_ab` at `(a, b)`.
    ///
    /// # Complexity
    /// - Time: amortized O(1)
    /// - Space: amortized O(1)
    ///
    /// # Panics
    /// Panics if `i >= n`, `j >= n`, `i == j`, `x00` or `x11` is `T::max_value()`, or
    /// `x00 + x11 > x01 + x10`.
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

    /// The minimum of the function, and a minimizer, with `true` for `x_i = 1`.
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

    /// Adds the edge `u -> v` of capacity `c` if `c > 0`.
    ///
    /// # Complexity
    /// - Time: amortized O(1)
    /// - Space: amortized O(1)
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
