/// The strongly connected components of a directed graph, in topological order.
///
/// # Definition
/// `G` is the directed graph on `[0, n)` with an edge `u -> v` for each `v` in `adjacency[u]`, and
/// `m` is its number of edges. `u` and `v` are strongly connected if each is reachable from the
/// other, and the classes of this equivalence relation are the components; `c` is their number.
/// Contracting each component to a vertex gives the condensation of `G`. The components are
/// numbered `0, ..., c - 1` in a topological order of the condensation: for every edge `u -> v`,
/// `component(u) <= component(v)`.
///
/// # Invariants
/// - `component[v]` is the number of the component of `v`.
/// - `members[start[k]..start[k + 1]]` is the component `k` in increasing order.
/// - `condensation[k]` is the `l != k` with an edge from the component `k` to the component `l`,
///   without duplicates.
///
/// # Complexity
/// - Space: O(n + m)
pub struct Scc {
    component: Box<[usize]>,
    start: Box<[usize]>,
    members: Box<[usize]>,
    condensation: Box<[Vec<usize>]>,
}

impl Scc {
    /// The strongly connected components of `G`.
    ///
    /// # Complexity
    /// - Time: O(n + m)
    /// - Space: O(n + m)
    ///
    /// # Panics
    /// Panics if some `v` in `adjacency[u]` satisfies `v >= n`.
    pub fn from_adjacency(adjacency: &[Vec<usize>]) -> Self {
        let n = adjacency.len();
        let mut order = vec![!0; n];
        let mut low = vec![0; n];
        let mut component = vec![!0; n];
        let mut visited = 0;
        let mut count = 0;
        let mut open = Vec::with_capacity(n);
        let mut stack = Vec::new();
        for root in 0..n {
            if order[root] != !0 {
                continue;
            }
            order[root] = visited;
            low[root] = visited;
            visited += 1;
            open.push(root);
            stack.push((root, 0));
            while let Some(&mut (v, ref mut i)) = stack.last_mut() {
                if *i < adjacency[v].len() {
                    let w = adjacency[v][*i];
                    *i += 1;
                    assert!(w < n, "vertex out of bounds: v={w}, n={n}");
                    if order[w] == !0 {
                        order[w] = visited;
                        low[w] = visited;
                        visited += 1;
                        open.push(w);
                        stack.push((w, 0));
                    } else if order[w] < n {
                        low[v] = low[v].min(order[w]);
                    }
                    continue;
                }
                stack.pop();
                if let Some(&(p, _)) = stack.last() {
                    low[p] = low[p].min(low[v]);
                }
                if low[v] == order[v] {
                    while let Some(u) = open.pop() {
                        order[u] = n;
                        component[u] = count;
                        if u == v {
                            break;
                        }
                    }
                    count += 1;
                }
            }
        }
        for c in &mut component {
            *c = count - 1 - *c;
        }
        let mut start = vec![0; count + 1];
        for &k in &component {
            start[k + 1] += 1;
        }
        for k in 0..count {
            start[k + 1] += start[k];
        }
        let mut members = vec![0; n];
        let mut next = start.clone();
        for v in 0..n {
            members[next[component[v]]] = v;
            next[component[v]] += 1;
        }
        let mut condensation = vec![Vec::new(); count];
        let mut mark = vec![!0; count];
        for k in 0..count {
            for &u in &members[start[k]..start[k + 1]] {
                for &v in &adjacency[u] {
                    let l = component[v];
                    if l != k && mark[l] != k {
                        mark[l] = k;
                        condensation[k].push(l);
                    }
                }
            }
        }
        Self {
            component: component.into(),
            start: start.into(),
            members: members.into(),
            condensation: condensation.into(),
        }
    }

    /// The number `component(v)` of the component of `v`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `v >= n`.
    pub fn component(&self, v: usize) -> usize {
        assert!(
            v < self.component.len(),
            "vertex out of bounds: v={v}, len={}",
            self.component.len()
        );
        self.component[v]
    }

    /// The component `k`, as its vertices in increasing order.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `k >= c`.
    pub fn group(&self, k: usize) -> &[usize] {
        assert!(
            k < self.len(),
            "component out of bounds: k={k}, len={}",
            self.len()
        );
        &self.members[self.start[k]..self.start[k + 1]]
    }

    /// The condensation of `G`, as a simple directed acyclic graph on `[0, c)`.
    ///
    /// # Definition
    /// The `k`-th list has the `l != k` with an edge from the component `k` to the component `l`,
    /// without duplicates, in an unspecified order.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn condensation(&self) -> &[Vec<usize>] {
        &self.condensation
    }

    /// The number `c` of components.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.start.len() - 1
    }

    /// Whether `n = 0`, that is `c = 0`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        self.component.is_empty()
    }
}
