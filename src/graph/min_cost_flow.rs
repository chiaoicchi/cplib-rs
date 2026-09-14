use crate::algebra::{Bounded, Zero};
use crate::graph::topological_sort::topological_sort;

/// A flow network with costs and a minimum cost flow algorithm (successive shortest paths).
///
/// # Definition
/// A flow network whose edges carry a capacity in `Cap` and a cost per unit of flow in `Cost`. The
/// cost of a flow is the sum over the edges of the flow times the cost. For an amount `x`, the
/// minimum cost of a flow of value `x` from `s` to `t` is a convex piecewise linear function
/// `g(x)`; the structure stores a flow, initially zero, and `slope` increases it along `g`, whose
/// breakpoints it returns.
///
/// The residual graph has, for each edge with capacity `c`, flow `f`, and cost `w`, a forward edge
/// of residual capacity `c - f` and cost `w`, and a reverse edge of residual capacity `f` and cost
/// `-w`. A flow is of minimum cost among the flows of its value if and only if the residual graph
/// has no cycle of negative cost; pushing along a shortest path of the residual graph keeps this
/// property, so pushing along shortest paths in turn traces `g`.
///
/// Costs may be negative, but the graph must have no cycle of negative cost. Shortest paths are
/// computed by Dijkstra's algorithm on the reduced costs `w + h(u) - h(v)`, where the potential `h`
/// keeps them non-negative; it is recomputed when needed.
///
/// # Invariants
/// - Edge `i` is stored as the pair of internal edges `2i` (forward) and `2i + 1` (reverse), with
///   `to[2i + 1]` its tail and `to[2i]` its head; `cost[2i + 1] = -cost[2i]`.
/// - `cap[e]` is the residual capacity of internal edge `e`; `cap[2i] + cap[2i + 1]` is the
///   capacity of edge `i`, and `cap[2i + 1]` is its flow.
/// - `adjacency[v]` lists the internal edges leaving `v`.
/// - After `slope`, `cost[e] + potential[from] - potential[to] >= 0` for every internal edge `e`
///   with `cap[e] > 0`.
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
    /// Constructs a network on `[0, n)` with no edges.
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

    /// Constructs a network on `[0, n)` with the given edges, so that edge `i` is `edges[i]`.
    /// Equivalent to `new(n)` followed by `add_edge` for each element in order.
    ///
    /// # Complexity
    /// - Time: O(n + m)
    /// - Space: O(n + m)
    ///
    /// # Panics
    /// Panics if some `from`, `to` is not less than `n`, or some capacity is negative.
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
    /// - Time: O(1) amortized
    /// - Space: O(1)
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

    /// Increases the stored flow from `s` to `t` to a maximum at minimum cost, and returns the
    /// increase of the value and of the cost.
    ///
    /// # Complexity
    /// - Time: O(F (n + m) log n), where `F` is the number of augmentations, plus the
    /// recomputation of the potential if needed.
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `s >= n`, `t >= n`, `s == t`, or the graph has a cycle of negative cost.
    pub fn flow(&mut self, s: usize, t: usize) -> (Cap, Cost) {
        self.flow_limit(s, t, Cap::max_value())
    }

    /// Increases the stored flow from `s` to `t` by at most `limit` at minimum cost, as much as
    /// possible, and returns the increase of the value and of the cost.
    ///
    /// # Complexity
    /// - Time: O(F (n + m) log n), plus the recomputation of the potential if needed
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

    /// Increases the stored flow from `s` to `t` by exactly `amount` at minimum cost, and returns
    /// the increase of the cost, or `None` if `amount` cannot be pushed; in that case the flow is
    /// still increased as much as possible.
    ///
    /// # Complexity
    /// - Time: O(F (n + m) log n), plus the recomputation of the potential if needed.
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `s >= n`, `t >= n`, `s == t`, or the graph has a cycle of negative cost.
    pub fn flow_exact(&mut self, s: usize, t: usize, amount: Cap) -> Option<Cost> {
        let (value, cost) = self.flow_limit(s, t, amount);
        (value == amount).then_some(cost)
    }

    /// Increases the stored flow from `s` to `t` to a maximum at minimum cost, and returns the
    /// breakpoints of `g` from `(0, 0)` to the end, relative to the flow before the call.
    ///
    /// # Complexity
    /// - Time: O(F (n + m) log n), plus the recomputation of the potential if needed
    /// - Space: O(n + F)
    ///
    /// # Panics
    /// Panics if `s >= n`, `t >= n`, `s == t`, or the graph has a cycle of negative cost.
    pub fn slope(&mut self, s: usize, t: usize) -> Vec<(Cap, Cost)> {
        self.slope_limit(s, t, Cap::max_value())
    }

    /// Increases the stored flow from `s` to `t` by at most `limit` at minimum cost, as much as
    /// possible, and returns the breakpoints of `g` from `(0, 0)` to the end, relative to the flow
    /// before the call. Consecutive segments of equal slope are merged.
    ///
    /// # Complexity
    /// - Time: O(F (n + m) log n), where `F` is the number of augmentations
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

    /// Recomputes the potential if some edge of positive residual capacity has negative reduced
    /// cost, so that all reduced costs are non-negative afterwards.
    ///
    /// # Complexity
    /// - Time: O(n + m) when residual graph is acyclic, O(nm) otherwise.
    /// - Space: O(n)
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
    /// # Complexity
    /// - Time: O((n + m) log n)
    /// - Space: O(n)
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

    /// Returns `(from, to, cap, flow, cost)` of edge `i`.
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

    /// Returns `n`, the number of vertices.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn num_vertices(&self) -> usize {
        self.adjacency.len()
    }

    /// Returns `m`, the number of edges.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn num_edges(&self) -> usize {
        self.to.len() >> 1
    }
}
