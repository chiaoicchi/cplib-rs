use crate::algebra::{Bounded, Zero};

/// A flow network with a maximum flow algorithm (Dinic).
///
/// # Definition
/// A flow network is a directed graph on `[0, n)` with `m` edges, whose edge `i` carries a capacity
/// `c(i)` in `Cap`. The structure stores a value `f(i)` in `[0, c(i)]` on each edge, initially `0`.
/// A flow from `s` to `t` is such an `f` with the total inflow equal to the total outflow at every
/// vertex other than `s` and `t`, and its value is the net outflows of `s`.
///
/// The residual graph has, for each edge `i` from `u` to `v`, a forward edge `u -> v` of residual
/// capacity `c(i) - f(i)` and a reverse edge `v -> u` of residual capacity `f(i)`.
///
/// # Invariants
/// - Edge `i` is stored as the pair of internal edges `2i` (forward) and `2i + 1` (reverse), with
///   `to[2i + 1]` its tail and `to[2i]` its head.
/// - `cap[e]` is the residual capacity of internal edge `e`, so that `cap[2i] + cap[2i + 1] = c(i)`
///   and `cap[2i + 1] = f(i)`.
/// - `adjacency[v]` lists the internal edges leaving `v`.
///
/// # Complexity
/// - Space: O(n + m)
pub struct MaxFlow<Cap> {
    to: Vec<usize>,
    cap: Vec<Cap>,
    adjacency: Vec<Vec<usize>>,
}

impl<Cap> MaxFlow<Cap> {
    /// The network on `[0, n)` with no edges.
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

impl<Cap: Copy + Ord + Zero + Bounded + std::ops::Add<Output = Cap> + std::ops::Sub<Output = Cap>>
    MaxFlow<Cap>
{
    /// The network on `[0, n)` whose edge `i` is `edges[i]`, as by `new(n)` followed for each edge
    /// in order.
    ///
    /// # Complexity
    /// - Time: O(n + m)
    /// - Space: O(n + m)
    ///
    /// # Panics
    /// Panics if some edge has `from >= n`, `to >= n` or `cap < Cap::zero()`.
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
    /// - Time: amortized O(1)
    /// - Space: amortized O(1)
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

    /// Pushes along paths from `s` to `t` in the residual graph until `t` in unreachable, and
    /// returns the total amount, as by `flow_limit(s, t, Cap::max_value())`.
    ///
    /// # Definition
    /// The push keeps the conservation at every vertex other than `s` and `t`, and raises the net
    /// outflow of `s` by the returned amount. If `f` was a flow from `s` to `t`, it becomes a
    /// maximum one.
    ///
    /// # Complexity
    /// - Time: O(n^2 m), and O(m min(√m, n^{2/3})) if every capacity is at most `1`
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `s >= n`, `t >= n`, or `s == t`.
    pub fn flow(&mut self, s: usize, t: usize) -> Cap {
        self.flow_limit(s, t, Cap::max_value())
    }

    /// Pushes along paths from `s` to `t` in the residual graph, as much as possible up to `limit`,
    /// and returns the total amount.
    ///
    /// # Definition
    /// The push keeps the conservation at every vertex other than `s` and `t`, and raises the net
    /// outflow of `s` by the returned amount, which is `limit` or the most that can be pushed.
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
        let mut queue = Vec::with_capacity(n);
        let mut path = Vec::new();
        while total < limit {
            dist.fill(!0);
            dist[s] = 0;
            queue.clear();
            queue.push(s);
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

    /// Pushes along one path from `s` to `t` of the level graph given by `dist`, found from `iter`,
    /// at most `limit`, and returns the amount, or `Cap::zero()` if no such path remains.
    ///
    /// # Complexity
    /// - Time: O(n + k), where `k` is the number of edges discarded from `iter`, which total `O(m)`
    ///   over a phase
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

    /// The edge `i`, as `(from, to, c(i), f(i))`.
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

    /// Whether each vertex is reachable from `s` in the residual graph.
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
}
