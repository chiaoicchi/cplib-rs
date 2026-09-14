/// Returns a topological order of a directed graph, and `None` if it has a cycle.
///
/// # Definition
/// A topological order of a directed graph `G` on `[0, n)` is a permutation `p` of
/// `[0, n)` such that `u` appears before `v` in `p` for every edge `u -> v`. It exists if and only
/// if `G` is acyclic. When several exist, which one is returned is unspecified.
///
/// # Contract
/// `adjacency[u]` is the list of `v` for the edges `u -> v` of `G`, with `v < n`.
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
    let mut stack: Vec<usize> = (0..n).filter(|i| degree[*i] == 0).collect();
    let mut order = Vec::with_capacity(n);
    while let Some(u) = stack.pop() {
        order.push(u);
        for &v in &adjacency[u] {
            degree[v] -= 1;
            if degree[v] == 0 {
                stack.push(v);
            }
        }
    }
    (order.len() == n).then_some(order)
}
