use crate::algebra::{Bounded, Zero};

/// A flow network with a maximum flow algorithm (Dinic).
///
/// # Definition
/// A flow network is a directed graph on `[0, n)` whose edges carry capacities in `Cap`. A flow
/// from `s` to `t` assigns to each edge a value in `[0, capacity]` such that, at every vertex other
/// than `s` and `t`, the total inflow equals the total outflow; its value is the total outflow of
/// `s`. The structure stores a flow, initially zero, and `flow` increases it to a maximum.
///
/// The residual graph has, for each edge `i` from `u` to `v` with capacity `c` and flow `f`, a
/// forward edge `u -> v` of residual capacity `c - f` and a reverse edge `v -> u` of residual
/// capacity `f`. Pushing along a path of the residual graph increases the flow; the flow is maximum
/// if and only if `t` is not reachable from `s` in the residual graph.
///
/// # Invariants
/// - Edge `i` is stored as the pair of internal edges `2i` (forward) and `2i + 1` (reverse), with
///   `to[2i + 1]` its tail and `to[2i]` its head.
/// - `cap[e]` is the residual capacity of internal edge `e`; `cap[2i] + cap[2i + 1]` is the
///   capacity of edge `i`, and `cap[2i + 1]` is its flow.
/// - `adjacency[v]` lists the internal edges leaving `v`.
///
/// # Complexity
/// - Space: O(n + m)
pub struct MaxFlow<Cap> {
    to: Vec<usize>,
    cap: Vec<Cap>,
    adjacency: Vec<Vec<usize>>,
}

impl<Cap: Copy + Ord + Zero + Bounded + std::ops::Add<Output = Cap> + std::ops::Sub<Output = Cap>>
    MaxFlow<Cap>
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
            adjacency: vec![Vec::new(); n],
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
    pub fn from_edges(n: usize, edges: &[(usize, usize, Cap)]) -> Self {
        let mut network = Self::new(n);
        for &(from, to, cap) in edges {
            network.add_edge(from, to, cap);
        }
        network
    }

    /// Adds an edge from `from` to `to` with capacity `cap`, and returns its index.
    ///
    /// # Complexity
    /// - Time: O(1) amortized
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `from >= n`, `to >= n`, or `cap < zero`.
    pub fn add_edge(&mut self, from: usize, to: usize, cap: Cap) -> usize {
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
        e / 2
    }

    /// Increases the stored flow from `s` to `t` to a maximum, and returns the increase.
    /// Equivalent to `flow_limit(s, t, Cap::max_value())`.
    ///
    /// # Complexity
    /// - Time: O(n^2 m); O(m sqrt(n)) when all capacities are `one`.
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `s >= n`, `t >= n`, or `s == t`.
    pub fn flow(&mut self, s: usize, t: usize) -> Cap {
        self.flow_limit(s, t, Cap::max_value())
    }

    /// Increases the stored flow from `s` to `t` by at most `limit`, as much as possible, and
    /// returns the increase.
    ///
    /// # Complexity
    /// - Time: O(n^2 m)
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `s >= n`, `t >= n`, or `s == t`.
    pub fn flow_limit(&mut self, s: usize, t: usize, limit: Cap) -> Cap {
        let n = self.adjacency.len();
        assert!(s < n, "vertex out of bounds: s={s}, n={n}");
        assert!(t < n, "vertex out of bounds: t={t}, n={n}");
        assert!(s != t, "source and sink must differ: s={s}");
        let mut total = Cap::zero();
        let mut dist = vec![!0usize; n];
        let mut iter = vec![0; n];
        let mut path = Vec::new();
        while total < limit {
            dist.fill(!0);
            dist[s] = 0;
            let mut queue = vec![s];
            let mut head = 0;
            while head < queue.len() && dist[t] == !0 {
                let v = queue[head];
                head += 1;
                for &e in &self.adjacency[v] {
                    let w = self.to[e];
                    if self.cap[e] > Cap::zero() && dist[w] == !0 {
                        dist[w] = dist[v] + 1;
                        queue.push(w);
                    }
                }
            }
            if dist[t] == !0 {
                break;
            }
            iter.fill(0);
            while total < limit {
                let pushed = self.augment(s, t, limit - total, &dist, &mut iter, &mut path);
                if pushed == Cap::zero() {
                    break;
                }
                total = total + pushed;
            }
        }
        total
    }

    /// Pushes along one `s`-`t` path of the level graph, from `iter`, and returns the amount;
    /// returns `zero` if no such path remains.
    ///
    /// # Complexity
    /// - Time: O(n + k), where `k` is the number of edges discarded from `iter`; amortized over a
    ///   phase, the discarded edges total `O(m)`.
    /// - Space: O(1)
    fn augment(
        &mut self,
        s: usize,
        t: usize,
        limit: Cap,
        dist: &[usize],
        iter: &mut [usize],
        path: &mut Vec<usize>,
    ) -> Cap {
        path.clear();
        let mut v = s;
        while v != t {
            let mut next = None;
            while iter[v] < self.adjacency[v].len() {
                let e = self.adjacency[v][iter[v]];
                let w = self.to[e];
                if self.cap[e] > Cap::zero() && dist[w] == dist[v] + 1 {
                    next = Some(e);
                    break;
                }
                iter[v] += 1;
            }
            match next {
                Some(e) => {
                    path.push(e);
                    v = self.to[e];
                }
                None => {
                    let Some(e) = path.pop() else {
                        return Cap::zero();
                    };
                    v = self.to[e ^ 1];
                    iter[v] += 1;
                }
            }
        }
        let mut d = limit;
        for &e in path.iter() {
            d = d.min(self.cap[e]);
        }
        for &e in path.iter() {
            self.cap[e] = self.cap[e] - d;
            self.cap[e ^ 1] = self.cap[e ^ 1] + d;
        }
        d
    }

    /// Returns `(from, to, cap, flow)` of edge `i`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i >= m`.
    pub fn get_edge(&self, i: usize) -> (usize, usize, Cap, Cap) {
        let m = self.num_edges();
        assert!(i < m, "edge out of bounds: i={i}, m={m}");
        let e = i << 1;
        (
            self.to[e ^ 1],
            self.to[e],
            self.cap[e] + self.cap[e ^ 1],
            self.cap[e ^ 1],
        )
    }

    /// Returns which vertices are reachable from `s` in the residual graph.
    ///
    /// Immediately after `flow(s, t)`, this is the `s` side of a minimum `s` - `t` cut: the edges
    /// from the reachable side to the other side are saturated, and their capacities sum to the
    /// value of the flow.
    ///
    /// # Complexity
    /// - Time: O(n + m)
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `s >= n`.
    pub fn reachable(&self, s: usize) -> Vec<bool> {
        let n = self.adjacency.len();
        assert!(s < n, "vertex out of bounds: s={s}, n={n}");
        let mut reached = vec![false; n];
        reached[s] = true;
        let mut stack = vec![s];
        while let Some(v) = stack.pop() {
            for &e in &self.adjacency[v] {
                let w = self.to[e];
                if self.cap[e] > Cap::zero() && !reached[w] {
                    reached[w] = true;
                    stack.push(w);
                }
            }
        }
        reached
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
