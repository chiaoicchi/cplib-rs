use crate::collections::dsu::Dsu;

/// A minimum spanning forest of an undirected weighted graph, as the edges it contains.
///
/// # Definition
/// `G` is the undirected multigraph on `[0, n)` whose `i`-th edge joins `u` and `v` with weight `w`
/// for `edges[i] = (u, v, w)`, and `m` is the length of `edges`. A spanning forest of `G` is a set
/// `F` of edges of `G` without cycles in which two vertices are connected if and only if they are
/// connected in `G`; when `G` is connected, it is a spanning tree. `F` is minimum if, for every
/// spanning forest `F'` and every `k`, the `k`-th least weight of `F` is at most that of `F'`.
/// Returns `in_forest` of length `m`, where `in_forest[i]` is whether the `i`-th edge is in `F` for
/// a minimum spanning forest `F`; which one is unspecified when several exist.
///
/// # Complexity
/// - Time: O(n + m log m)
/// - Space: O(n + m)
///
/// # Panics
/// Panics if some `(u, v, w)` in `edges` satisfies `u >= n` or `v >= n`.
pub fn minimum_spanning_tree<W: Ord>(n: usize, edges: &[(usize, usize, W)]) -> Vec<bool> {
    let mut order: Vec<usize> = (0..edges.len()).collect();
    order.sort_unstable_by(|&i, &j| edges[i].2.cmp(&edges[j].2));
    let mut dsu = Dsu::new(n);
    let mut in_forest = vec![false; edges.len()];
    for i in order {
        let (u, v, _) = edges[i];
        assert!(u < n, "vertex out of bounds: u={u}, n={n}");
        assert!(v < n, "vertex out of bounds: v={v}, n={n}");
        in_forest[i] = dsu.unite(u, v);
    }
    in_forest
}
