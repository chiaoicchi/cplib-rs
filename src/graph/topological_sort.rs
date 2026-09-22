/// A topological order of a directed graph `G`, or `None` if `G` has a cycle.
///
/// # Definition
/// `G` is the directed graph on `[0, n)` with an edge `u -> v` for each `v` in `adjacency[u]`, and
/// `m` is its number of edges. A topological order of `G` is a permutation `p` of `[0, n)` in
/// which `u` comes before `v` for every edge `u -> v`. When several exist, which one is returned is
/// unspecified.
///
/// # Complexity
/// - Time: O(n + m)
/// - Space: O(n)
///
/// # Panics
/// Panics if some `v` in `adjacency[u]` satisfies `v >= n`.
pub fn topological_sort(adjacency: &[Vec<usize>]) -> Option<Vec<usize>> {
    let n = adjacency.len();
    let mut degree = vec![0usize; n];
    for edges in adjacency {
        for &v in edges {
            assert!(v < n, "vertex out of bounds: v={v}, n={n}");
            degree[v] += 1;
        }
    }
    let mut order = Vec::with_capacity(n);
    order.extend((0..n).filter(|&i| degree[i] == 0));
    let mut i = 0;
    while i < order.len() {
        let u = order[i];
        i += 1;
        for &v in &adjacency[u] {
            degree[v] -= 1;
            if degree[v] == 0 {
                order.push(v);
            }
        }
    }
    (order.len() == n).then_some(order)
}
