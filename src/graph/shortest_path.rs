use std::borrow::Borrow;

use crate::algebra::Zero;

/// A shortest-path forest of a directed graph from a set of sources.
///
/// # Definition
/// `G` is a directed graph on `[0, n)` with weighted edges and `S` is a set of sources, both given
/// to the constructor. For a vertex `v`, `d(v)` is the minimum total weight of a walk in `G` from a
/// vertex of `S` to `v`, and is undefined if no such walk exists. `p(v)` is a vertex such that `G`
/// has an edge from `p(v)` to `v` of weight `w` with `d(p(v)) + w = d(v)`, and is undefined if
/// `d(v)` is undefined or if `v` is in `S` and `d(v) = 0`. Following `p` from a vertex `v` with
/// `d(v)` defined ends at a vertex of `S`. When several choices of `p` exist, which one is taken is
/// unspecified.
///
/// # Invariants
/// - `dist[v]` is `d(v)`, or `None` if undefined.
/// - `prev[v]` is `p(v)`, or `None` if undefined.
///
/// # Complexity
/// - Space: O(n)
pub struct ShortestPathTree<W> {
    dist: Vec<Option<W>>,
    prev: Vec<Option<usize>>,
}

impl ShortestPathTree<usize> {
    /// The shortest-path forest of an unweighted graph, by breadth-first search.
    ///
    /// # Definition
    /// `G` has an edge of weight `1` from `v` to each `u` in `edges(v)`, and `S` is the set of
    /// `sources`. `m` is the total nember of items of `edges(v)` over all `v`. `edges(v)` is called
    /// at most once for each `v`.
    ///
    /// # Complexity
    /// - Time: O(n + m)
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if some `s` in `sources` or some `u` in some `edges(v)` satisfies `s >= n` or
    /// `u >= n`.
    pub fn bfs<I>(n: usize, sources: &[usize], mut edges: impl FnMut(usize) -> I) -> Self
    where
        I: IntoIterator,
        I::Item: Borrow<usize>,
    {
        let mut dist = vec![None; n];
        let mut prev = vec![None; n];
        let mut queue = std::collections::VecDeque::new();
        for &s in sources {
            assert!(s < n, "vertex out of bounds: s={s}, n={n}");
            if dist[s].is_none() {
                dist[s] = Some(0);
                queue.push_back(s);
            }
        }

        while let Some(v) = queue.pop_front() {
            let d = dist[v].unwrap() + 1;
            for u in edges(v) {
                let u = *u.borrow();
                assert!(u < n, "vertex out of bounds: u={u}, n={n}");
                if dist[u].is_none() {
                    dist[u] = Some(d);
                    prev[u] = Some(v);
                    queue.push_back(u);
                }
            }
        }
        Self { dist, prev }
    }

    /// The shortest-path forest of a graph with weights `0` and `1`, by 0-1 breadth-first search.
    ///
    /// # Definition
    /// `G` has an edge of weight `w` from `v` to `u` for each `(u, w)` in `edges(v)`, and `S` is
    /// the set of `sources`. `m` is the total number of items of `edges(v)` over all `v`.
    /// `edges(v)` is called at most once for each `v`.
    ///
    /// # Complexity
    /// - Time: O(n + m)
    /// - Space: O(n + m)
    ///
    /// # Panics
    /// Panics if some `s` in `sources` or some `(u, w)` in some `edges(v)` satisfies `s >= n`,
    /// `u >= n` or `w > 1`.
    pub fn zero_one_bfs<I>(n: usize, sources: &[usize], mut edges: impl FnMut(usize) -> I) -> Self
    where
        I: IntoIterator,
        I::Item: Borrow<(usize, usize)>,
    {
        let mut dist = vec![None; n];
        let mut prev = vec![None; n];
        let mut deque = std::collections::VecDeque::new();
        for &s in sources {
            assert!(s < n, "vertex out of bounds: s={s}, n={n}");
            if dist[s].is_none() {
                dist[s] = Some(0);
                deque.push_back((0, s));
            }
        }

        while let Some((d, v)) = deque.pop_front() {
            if dist[v] != Some(d) {
                continue;
            }
            for e in edges(v) {
                let &(u, w) = e.borrow();
                assert!(u < n, "vertex out of bounds: u={u}, n={n}");
                assert!(w <= 1, "weight must be 0 or 1: w={w}");
                if dist[u].is_none_or(|x| d + w < x) {
                    dist[u] = Some(d + w);
                    prev[u] = Some(v);
                    if w == 0 {
                        deque.push_front((d, u));
                    } else {
                        deque.push_back((d + 1, u));
                    }
                }
            }
        }
        Self { dist, prev }
    }
}

impl<W: Copy + Ord + std::ops::Add<Output = W> + Zero> ShortestPathTree<W> {
    /// The shortest-path forest of a graph with non-negative weights, by Dijkstra's algorithm.
    ///
    /// # Definition
    /// `G` has an edge of weight `w` from `v` to `u` for each `(u, w)` in `edges(v)`, and `S` is
    /// the set of `sources`. `m` is the total number of items of `edges(v)` over all `v`.
    /// `edges(v)` is called at most once for each `v`.
    ///
    /// # Contract
    /// - `d(v) + w` does not overflow `W` for every edge `(u, w)` in `edges(v)`.
    ///
    /// # Complexity
    /// - Time: O(n + m log m)
    /// - Space: O(n + m)
    ///
    /// # Panics
    /// Panics if some `s` in `sources` or some `(u, w)` in some `edges(v)` satisfies `s >= n`,
    /// `u >= n` or `w < W::zero()`.
    pub fn dijkstra<I>(n: usize, sources: &[usize], mut edges: impl FnMut(usize) -> I) -> Self
    where
        I: IntoIterator,
        I::Item: std::borrow::Borrow<(usize, W)>,
    {
        let mut dist = vec![None; n];
        let mut prev = vec![None; n];
        let mut heap = std::collections::BinaryHeap::new();
        for &s in sources {
            assert!(s < n, "vertex out of bounds: s={s}, n={n}");
            if dist[s].is_none() {
                dist[s] = Some(W::zero());
                heap.push(std::cmp::Reverse((W::zero(), s)));
            }
        }

        while let Some(std::cmp::Reverse((d, v))) = heap.pop() {
            if dist[v] != Some(d) {
                continue;
            }
            for e in edges(v) {
                let &(u, w) = e.borrow();
                assert!(u < n, "vertex out of bounds: u={u}, n={n}");
                assert!(w >= W::zero(), "weight must be non-negative");
                if dist[u].is_none_or(|x| d + w < x) {
                    dist[u] = Some(d + w);
                    prev[u] = Some(v);
                    heap.push(std::cmp::Reverse((d + w, u)));
                }
            }
        }
        Self { dist, prev }
    }

    /// The shortest-path forest of a graph with possibly negative weights, by the Bellman-Ford
    /// algorithm, or `None` if a cycle of negative total weight is reachable from the sources.
    ///
    /// # Definition
    /// `G` has an edge of weight `w` from `v` to `u` for each `(u, w)` in `edges(v)`, and `S` is
    /// the set of `sources`. `m` is the total number of items of `edges(v)` over all `v`. Returns
    /// `None` if and only if `G` has a cycle of negative total weight reachable from a vertex of
    /// `S`.
    ///
    /// # Contract
    /// - `edges(v)` yields the same items every time it is called.
    /// - No sum of weights along a walk from a vertex of `S` with at most `n` edges overflows `W`.
    ///
    /// # Complexity
    /// - Time: O(n(n + m))
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if some `s` in `sources` or some `(u, w)` in some `edges(v)` satisfies `s >= n` or
    /// `u >= n`.
    pub fn bellman_ford<I>(
        n: usize,
        sources: &[usize],
        mut edges: impl FnMut(usize) -> I,
    ) -> Option<Self>
    where
        I: IntoIterator,
        I::Item: std::borrow::Borrow<(usize, W)>,
    {
        let mut dist = vec![None; n];
        let mut prev = vec![None; n];
        for &s in sources {
            assert!(s < n, "vertex out of bounds: s={s}, n={n}");
            dist[s] = Some(W::zero());
        }
        for _ in 0..=n {
            let mut updated = false;
            for v in 0..n {
                let Some(d) = dist[v] else {
                    continue;
                };
                for e in edges(v) {
                    let &(u, w) = e.borrow();
                    assert!(u < n, "vertex out of bounds: u={u}, n={n}");
                    if dist[u].is_none_or(|x| d + w < x) {
                        dist[u] = Some(d + w);
                        prev[u] = Some(v);
                        updated = true;
                    }
                }
            }
            if !updated {
                return Some(Self { dist, prev });
            }
        }
        None
    }
}

impl<W: Copy> ShortestPathTree<W> {
    /// `d(v)`, or `None` if `d(v)` is undefined.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `v >= n`.
    pub fn dist(&self, v: usize) -> Option<W> {
        self.dist[v]
    }

    /// `p(v)`, or `None` if `p(v)` is undefined.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `v >= n`.
    pub fn prev(&self, v: usize) -> Option<usize> {
        self.prev[v]
    }

    /// The vertices `v_0, ..., v_k` of a shortest walk from a vertex of `S` to `t`, or `None` if
    /// `d(t)` is undefined.
    ///
    /// # Definition
    /// `v_k = t`, `v_{i-1} = p(v_i)` for `1 <= i <= k`, and `p(v_0)` is undefined.
    ///
    /// # Complexity
    /// - Time: O(k)
    /// - Space: O(k)
    ///
    /// # Panics
    /// Panics if `t >= n`.
    pub fn path(&self, t: usize) -> Option<Vec<usize>> {
        self.dist[t]?;
        let mut path = vec![t];
        let mut v = t;
        while let Some(u) = self.prev[v] {
            path.push(u);
            v = u;
        }
        path.reverse();
        Some(path)
    }
}
