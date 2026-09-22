use crate::algebra::{Bounded, Zero};
use crate::graph::topological_sort::topological_sort;

/// A flow network with costs and a minimum cost flow algorithm (successive shortest paths).
///
/// # Definition
/// A flow network is a directed graph on `[0, n)` with `m` edges, whose edge `i` carries a capacity
/// `c(i)` in `Cap` and a cost `w(i)` per unit in `Cost`, which may be negative. The structure
/// stores a value `f(i)` in `[0, c(i)]` on each edge, initially `0`, and the cost of `f` is
/// `Σ_i f(i) w(i)`. A flow from `s` to `t` is such an `f` with the total inflow equal to the total
/// outflow at every vertex other than `s` and `t`, and its value is the net outflow of `s`. `g(x)`
/// is the minimum cost of a flow of value `x` from `s` to `t`.
///
/// The residual graph has, for each edge `i` from `u` to `v`, a forward edge `u -> v` of residual
/// capacity `c(i) - f(i)` and cost `w(i)`, and a reverse edge `v -> u` of residual capacity `f(i)`
/// and cost `-w(i)`. Only the edges of positive residual capacity are in it.
///
/// # Invariants
/// - Edge `i` is stored as the pair of internal edges `2i` (forward) and `2i + 1` (reverse), with
///   `to[2i + 1]` its tail and `to[2i]` its head, and `cost[2i] = -cost[2i + 1] = w(i)`.
/// - `cap[e]` is the residual capacity of internal edge `e`, so that `cap[2i] + cap[2i + 1] = c(i)`
///   and `cap[2i + 1] = f(i)`.
/// - `adjacency[v]` lists the internal edges leaving `v`.
/// - `cost[e] + potential[to[e^1]] - potential[to[e]] >= 0` for every internal edge `e` with
///   `cap[e] > 0`, except the edges added since the last push.
///
/// # Complexity
/// - Space: O(n + m)
pub struct MinCostFlow<Cap, Cost> {
    to: Vec<usize>,
    cap: Vec<Cap>,
    cost: Vec<Cost>,
    adjacency: Vec<Vec<usize>>,
    potential: Vec<Cost>,
}

impl<
    Cap: Copy + Ord + Zero + Bounded + std::ops::Add<Output = Cap> + std::ops::Sub<Output = Cap>,
    Cost: Copy
        + Ord
        + Zero
        + Bounded
        + From<Cap>
        + std::ops::Add<Output = Cost>
        + std::ops::Sub<Output = Cost>
        + std::ops::Mul<Output = Cost>,
> MinCostFlow<Cap, Cost>
{
    /// The network on `[0, n)` with no edges.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    pub fn new(n: usize) -> Self {
        Self {
            to: Vec::new(),
            cap: Vec::new(),
            cost: Vec::new(),
            adjacency: vec![Vec::new(); n],
            potential: vec![Cost::zero(); n],
        }
    }

    /// The network on `[0, n)` whose edge `i` is `edges[i]`, as by `new(n)` followrd by `add_edge`
    /// for each edge in order.
    ///
    /// # Complexity
    /// - Time: O(n + m)
    /// - Space: O(n + m)
    ///
    /// # Panics
    /// Panics if some edge has `from >= n`, `to >= n` or `cap < Cap::zero()`.
    pub fn from_edges(n: usize, edges: &[(usize, usize, Cap, Cost)]) -> Self {
        let mut network = Self::new(n);
        for &(from, to, cap, cost) in edges {
            network.add_edge(from, to, cap, cost);
        }
        network
    }

    /// Adds an edge from `from` to `to` with capacity `cap` and cost `cost`, and returns its index.
    ///
    /// # Complexity
    /// - Time: amortized O(1)
    /// - Space: amortized O(1)
    ///
    /// # Panics
    /// Panics if `from >= n`, `to >= n` or `cap < zero`.
    pub fn add_edge(&mut self, from: usize, to: usize, cap: Cap, cost: Cost) -> usize {
        let n = self.adjacency.len();
        assert!(from < n, "vertex out of bounds: from={from}, n={n}");
        assert!(to < n, "vertex out of bounds: to={to}, n={n}");
        assert!(cap >= Cap::zero(), "capacity must be non-negative");
        let e = self.to.len();
        self.adjacency[from].push(e);
        self.adjacency[to].push(e ^ 1);
        self.to.push(to);
        self.to.push(from);
        self.cap.push(cap);
        self.cap.push(Cap::zero());
        self.cost.push(cost);
        self.cost.push(Cost::zero() - cost);
        e / 2
    }

    /// Pushes along shortest paths from `s` to `t` in the residual graph until `t` is unreachable,
    /// and returns the total amount and its cost, as by `flow_limit(s, t, Cap::max_value())`.
    ///
    /// # Definition
    /// The push keeps the conservation at every vertex other than `s` and `t`, and raises the net
    /// outflow of `s` by the returned amount. If `f` was a minimum cost flow from `s` to `t` of
    /// value `v`, it becomes one of the maximum value `v'`, and the returned cost is
    /// `g(v') - g(v)`.
    ///
    /// # Complexity
    /// - Time: O((F + 1)(n + m) log n), where `F` is the number of pushed along, plus O(nm) if an
    ///   edge of the residual graph has a negative reduced cost, as at the first call with a
    ///   negative cost
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `s >= n`, `t >= n`, `s == t`, or the graph has a cycle of negative cost.
    pub fn flow(&mut self, s: usize, t: usize) -> (Cap, Cost) {
        self.flow_limit(s, t, Cap::max_value())
    }

    /// Pushes along shortest paths from `s` to `t` in the residual graph, as much as possible up to
    /// `limit`, and returns the total amount and its cost.
    ///
    /// # Definition
    /// The push keeps the conservation at every vertex other than `s` and `t`, and raises the net
    /// outflow of `s` by the returned amount `x`. If `f` was a minimum cost flow from `s` to `t` of
    /// value `v`, it becomes one of value `v + x`, and the returned cost is `g(v + x)`.
    ///
    /// # Complexity
    /// - Time: O((F + 1)(n + m) log n), where `F` is the number of paths pushed along, plus O(nm)
    ///   if an edge of the residual graph has a negative cost
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `s >= n`, `t >= n`, `s == t`, or the graph has a cycle of negative cost.
    pub fn flow_limit(&mut self, s: usize, t: usize, limit: Cap) -> (Cap, Cost) {
        let n = self.adjacency.len();
        assert!(s < n, "vertex out of bounds: s={s}, n={n}");
        assert!(t < n, "vertex out of bounds: t={t}, n={n}");
        assert!(s != t, "source and sink must differ: s={s}");
        self.restore_potential();
        let mut value = Cap::zero();
        let mut total = Cost::zero();
        let mut dist = vec![Cost::max_value(); n];
        let mut prev_edge = vec![!0; n];
        let mut done = vec![false; n];
        while value < limit {
            let Some((d, slope)) =
                self.augment(s, t, limit - value, &mut dist, &mut prev_edge, &mut done)
            else {
                break;
            };
            value = value + d;
            total = total + slope * Cost::from(d);
        }
        (value, total)
    }

    /// Pushes `amount` as by `flow_limit(s, t, amount)`, and returns its cost, or `None` if less
    /// than `amount` is pushed.
    ///
    /// # Complexity
    /// - Time: O((F + 1)(n + m) log n), where `F` is the number of paths pushed along, plus O(nm)
    ///   if an edge of the residual graph has a negative reduced cost
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `s >= n`, `t >= n`, `s == t`, or the graph has a cycle of negative cost.
    pub fn flow_exact(&mut self, s: usize, t: usize, amount: Cap) -> Option<Cost> {
        let (value, cost) = self.flow_limit(s, t, amount);
        (value == amount).then_some(cost)
    }

    /// Pushes as by `flow(s, t)` and returns the breakpoints of the cost against the amount, as by
    /// `slope_limit(s, t, Cap::max_value())`.
    ///
    /// # Complexity
    /// - Time: O((F + 1)(n + m) log n), where `F` is the number of paths pushed along, plus O(nm)
    ///   if an edge of the residual graph has a negative reduced cost
    /// - Space: O(n + F)
    ///
    /// # Panics
    /// Panics if `s >= n`, `t >= n`, `s == t`, or the graph has a cycle of negative cost.
    pub fn slope(&mut self, s: usize, t: usize) -> Vec<(Cap, Cost)> {
        self.slope_limit(s, t, Cap::max_value())
    }

    /// Pushes as by `flow_limit(s, t, limit)`, and returns the breakpoints of the cost against the
    /// amount.
    ///
    /// # Definition
    /// The breakpoints are the points `(x, y)`, from `(0, 0)` to the returned amount, where `y` is
    /// the cost of pushing `x` and the slope changes. If `f` was a minimum cost flow from `s` to `t`
    /// of value `v`, then `y = g(v + x) - g(v)`.
    ///
    /// # Complexity
    /// - Time: O((F + 1)(n + m) log n), where `F` is the number of paths pushed along, plus O(nm)
    ///   if an edge of the residual graph has a negative reduced cost
    /// - Space: O(n + F)
    ///
    /// # Panics
    /// Panics if `s >= n`, `t >= n`, `s == t`, or the graph has a cycle of negative cost.
    pub fn slope_limit(&mut self, s: usize, t: usize, limit: Cap) -> Vec<(Cap, Cost)> {
        let n = self.adjacency.len();
        assert!(s < n, "vertex out of bounds: s={s}, n={n}");
        assert!(t < n, "vertex out of bounds: t={t}, n={n}");
        assert!(s != t, "source and sink must differ: s={s}");
        self.restore_potential();
        let mut breakpoints = vec![(Cap::zero(), Cost::zero())];
        let mut value = Cap::zero();
        let mut total = Cost::zero();
        let mut prev_slope = None;
        let mut dist = vec![Cost::max_value(); n];
        let mut prev_edge = vec![!0; n];
        let mut done = vec![false; n];
        while value < limit {
            let Some((d, slope)) =
                self.augment(s, t, limit - value, &mut dist, &mut prev_edge, &mut done)
            else {
                break;
            };
            value = value + d;
            total = total + slope * Cost::from(d);
            if prev_slope == Some(slope) {
                breakpoints.pop();
            }
            breakpoints.push((value, total));
            prev_slope = Some(slope);
        }
        breakpoints
    }

    /// Recomputes the potential if an edge of the residual graph has a negative reduced
    /// cost, so that every reduced cost is non-negative.
    ///
    /// # Complexity
    /// - Time: O(m), plus O(n + m) for the recomputation if the residual graph is acyclic and
    ///   O(nm) otherwise
    /// - Space: O(n + m)
    fn restore_potential(&mut self) {
        let n = self.adjacency.len();
        let live = (0..self.to.len()).filter(|&e| self.cap[e] > Cap::zero());
        if live.clone().all(|e| {
            self.cost[e] + self.potential[self.to[e ^ 1]] - self.potential[self.to[e]]
                >= Cost::zero()
        }) {
            return;
        }
        let mut adjacency = vec![Vec::new(); n];
        for e in live.clone() {
            adjacency[self.to[e ^ 1]].push(self.to[e]);
        }

        let mut h = vec![Cost::zero(); n];
        if let Some(order) = topological_sort(&adjacency) {
            for &v in &order {
                for &e in &self.adjacency[v] {
                    if self.cap[e] > Cap::zero() {
                        let w = self.to[e];
                        h[w] = h[w].min(h[v] + self.cost[e]);
                    }
                }
            }
        } else {
            for round in 0..n {
                let mut relaxed = false;
                for e in live.clone() {
                    let (u, v) = (self.to[e ^ 1], self.to[e]);
                    if h[u] + self.cost[e] < h[v] {
                        h[v] = h[u] + self.cost[e];
                        relaxed = true;
                    }
                }
                if !relaxed {
                    break;
                }
                assert!(round + 1 < n, "the graph has a cycle of negative cost");
            }
        }
        self.potential = h;
    }

    /// Pushes along one shortest path from `s` to `t` in the residual graph, at most `limit`, and
    /// returns the amount and the cost per unit; `None` if `t` is unreachable.
    ///
    /// # Contract
    /// Every reduced cost of the residual graph is non-negative for `potential`.
    ///
    /// # Complexity
    /// - Time: O((n + m) log n)
    /// - Space: O(n + m)
    fn augment(
        &mut self,
        s: usize,
        t: usize,
        limit: Cap,
        dist: &mut [Cost],
        prev_edge: &mut [usize],
        done: &mut [bool],
    ) -> Option<(Cap, Cost)> {
        dist.fill(Cost::max_value());
        done.fill(false);
        dist[s] = Cost::zero();
        let mut heap = std::collections::BinaryHeap::new();
        heap.push(std::cmp::Reverse((Cost::zero(), s)));
        while let Some(std::cmp::Reverse((d, v))) = heap.pop() {
            if done[v] {
                continue;
            }
            done[v] = true;
            if v == t {
                break;
            }
            for &e in &self.adjacency[v] {
                let w = self.to[e];
                if self.cap[e] == Cap::zero() || done[w] {
                    continue;
                }
                let nd = d + self.cost[e] + self.potential[v] - self.potential[w];
                if nd < dist[w] {
                    dist[w] = nd;
                    prev_edge[w] = e;
                    heap.push(std::cmp::Reverse((nd, w)));
                }
            }
        }
        if !done[t] {
            return None;
        }
        for v in 0..self.adjacency.len() {
            if done[v] {
                self.potential[v] = self.potential[v] + dist[v] - dist[t];
            }
        }
        let slope = self.potential[t] - self.potential[s];
        let mut d = limit;
        let mut v = t;
        while v != s {
            let e = prev_edge[v];
            d = d.min(self.cap[e]);
            v = self.to[e ^ 1];
        }
        let mut v = t;
        while v != s {
            let e = prev_edge[v];
            self.cap[e] = self.cap[e] - d;
            self.cap[e ^ 1] = self.cap[e ^ 1] + d;
            v = self.to[e ^ 1];
        }
        Some((d, slope))
    }

    /// The edge `i`, as `(from, to, c(i), f(i), w(i))`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i >= m`.
    pub fn get_edge(&self, i: usize) -> (usize, usize, Cap, Cap, Cost) {
        assert!(
            i < self.num_edges(),
            "edge out of bounds: i={i}, m={}",
            self.num_edges()
        );
        let e = i << 1;
        (
            self.to[e ^ 1],
            self.to[e],
            self.cap[e] + self.cap[e ^ 1],
            self.cap[e ^ 1],
            self.cost[e],
        )
    }

    /// The number `n` of vertices.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn num_vertices(&self) -> usize {
        self.adjacency.len()
    }

    /// The number `m` of edges.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn num_edges(&self) -> usize {
        self.to.len() >> 1
    }
}
