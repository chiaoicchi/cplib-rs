/// The strongly connected components of a directed graph, in topological order.
///
/// # Definition
/// Two vertices `u`, `v` of a directed graph `G` on `[0, n)` are strongly connected if each is
/// reachable from the other. This is an equivalence relation, and its classes are the strongly
/// connected components. Contracting each component to a vertex gives the condensation of `G`,
/// which is a directed acyclic graph. Components are numbered `0, ..., c - 1` in a topological
/// order of the condensation: if there is an edge from `u` to `v` then
/// `component(u) <= component(v)`.
///
/// # Contract
/// `adjacency[u]` is the list of `v` for the edges `u -> v` of `G`, with `v < n`.
///
/// # Invariants
/// - `component[v] < c` for all `v`, and every value in `[0, c)` is taken.
/// - `groups[k]` is the list of vertices `v` with `component[v] == k`, in increasing order.
/// - For every edge `u -> v`, `component[u] <= component[v]`.
/// - `condensation[k]` lists distinct `l` with `l > k`.
///
/// # Complexity
/// - Space: O(n + m)
pub struct Scc {
    component: Box<[usize]>,
    groups: Box<[Vec<usize>]>,
    condensation: Box<[Vec<usize>]>,
}

impl Scc {
    /// Constructs the strongly connected components of `adjacency`.
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
        let mut groups = vec![Vec::new(); count];
        for v in 0..n {
            groups[component[v]].push(v);
        }
        let mut condensation = vec![Vec::new(); count];
        let mut mark = vec![!0; count];
        for (k, group) in groups.iter().enumerate() {
            for &u in group {
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
            groups: groups.into(),
            condensation: condensation.into(),
        }
    }

    /// Returns the component of `v`.
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

    /// Returns the vertices of component `k`, in increasing order.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `k >= c`.
    pub fn group(&self, k: usize) -> &[usize] {
        assert!(
            k < self.groups.len(),
            "component out of bounds: k={k}, len={}",
            self.groups.len()
        );
        &self.groups[k]
    }

    /// Returns the condensation of `adjacency` as a simple directed acyclic graph on `[0, c)`:
    /// `dag[k]` lists the `l` such that some edge goes from component `k` to component `l`,
    /// without duplicates, in the order first encountered Self-loops are omitted.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn condensation(&self) -> &[Vec<usize>] {
        &self.condensation
    }

    /// Returns `c`, the number of components.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.groups.len()
    }

    /// Returns `true` if the graph has no vertices.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        self.component.is_empty()
    }
}
