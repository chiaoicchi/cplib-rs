/// Calls `f` once for each triangle of an undirected graph.
///
/// # Definition
/// `G` is the simple undirected graph on `[0, n)` whose edges are the pairs `{u, v}` with `(u, v)`
/// in `edges`, and `m` is the length of `edges`. A triangle is a set `{a, b, c}` of three vertices
/// joined pairwise by edges of `G`. `f(a, b, c)` is called exactly once for each triangle; the
/// order of `a, b, c` within a call and the order of the calls are unspecified.
///
/// # Complexity
/// - Time: O(n + m√m)
/// - Space: O(n + m)
///
/// # Panics
/// Panics if some `(u, v)` in `edges` satisfies `u >= n`, `v >= n` or `u = v`.
pub fn for_each_triangle(
    n: usize,
    edges: &[(usize, usize)],
    mut f: impl FnMut(usize, usize, usize),
) {
    let mut simple: Vec<(usize, usize)> = edges
        .iter()
        .map(|&(u, v)| {
            assert!(u < n, "vertex out of bounds: u={u}, n={n}");
            assert!(v < n, "vertex out of bounds: v={v}, n={n}");
            assert!(u != v, "self loop: u={u}");
            (u.min(v), u.max(v))
        })
        .collect();
    simple.sort_unstable();
    simple.dedup();
    let mut degree = vec![0; n];
    for &(u, v) in &simple {
        degree[u] += 1;
        degree[v] += 1;
    }
    let mut out = vec![vec![]; n];
    for &(u, v) in &simple {
        if (degree[u], u) < (degree[v], v) {
            out[u].push(v);
        } else {
            out[v].push(u);
        }
    }
    let mut mark = vec![!0; n];
    for v in 0..n {
        for &w in &out[v] {
            mark[w] = v;
        }
        for &u in &out[v] {
            for &w in &out[u] {
                if mark[w] == v {
                    f(u, v, w);
                }
            }
        }
    }
}

/// The number of `C3` containing each edge of an undirected multigraph.
///
/// # Definition
/// `G` is the undirected multigraph on `[0, n)` whose `i`-th edge joins `u` and `v` for
/// `edges[i] = (u, v)`, and `m` is the length of `edges`. A `C3` of `G` is a set of three edges of
/// `G` forming a subgraph isomorphic to the cycle graph `C3`, that is, edges joining `{a, b}`,
/// `{b, c}` and `{c, a}` for three distinct vertices `a, b, c`. Returns `c` of length `m`, where
/// `c[i]` is the number of `C3` of `G` containing the `i`-th edge.
///
/// # Complexity
/// - Time: O(n + m√m)
/// - Space: O(n + m)
///
/// # Panics
/// Panics if some `(u, v)` in `edges` satisfies `u >= n`, `v >= n` or `u = v`.
pub fn count_c3(n: usize, edges: &[(usize, usize)]) -> Vec<usize> {
    let m = edges.len();
    let mut sorted: Vec<(usize, usize, usize)> = edges
        .iter()
        .enumerate()
        .map(|(i, &(u, v))| {
            assert!(u < n, "vertex out of bounds: u={u}, n={n}");
            assert!(v < n, "vertex out of bounds: v={v}, n={n}");
            assert!(u != v, "self loop: u={u}");
            (u.min(v), u.max(v), i)
        })
        .collect();
    sorted.sort_unstable();
    let mut class = vec![0; m];
    let mut ends = vec![];
    let mut mult = vec![];
    for (k, &(u, v, i)) in sorted.iter().enumerate() {
        if k == 0 || (sorted[k - 1].0, sorted[k - 1].1) != (u, v) {
            ends.push((u, v));
            mult.push(0);
        }
        class[i] = ends.len() - 1;
        *mult.last_mut().unwrap() += 1;
    }
    let mut degree = vec![0; n];
    for &(u, v) in &ends {
        degree[u] += 1;
        degree[v] += 1;
    }
    let mut out = vec![vec![]; n];
    for (e, &(u, v)) in ends.iter().enumerate() {
        if (degree[u], u) < (degree[v], v) {
            out[u].push((v, e));
        } else {
            out[v].push((u, e));
        }
    }
    let mut count = vec![0; ends.len()];
    let mut mark = vec![!0; n];
    let mut via = vec![0; n];
    for v in 0..n {
        for &(w, e) in &out[v] {
            mark[w] = v;
            via[w] = e;
        }
        for &(u, e1) in &out[v] {
            for &(w, e2) in &out[u] {
                if mark[w] == v {
                    let e3 = via[w];
                    count[e1] += mult[e2] * mult[e3];
                    count[e2] += mult[e3] * mult[e1];
                    count[e3] += mult[e1] * mult[e2];
                }
            }
        }
    }
    class.iter().map(|&e| count[e]).collect()
}

/// The number of `C4` containing each edge of an undirected multigraph.
///
/// # Definition
/// `G` is the undirected multigraph on `[0, n)` whose `i`-th edge joins `u` and `v` for
/// `edges[i] = (u, v)`, and `m` is the length of `edges`. A `C4` of `G` is a set of four edges of
/// `G` forming a subgraph isomorphic to the cycle graph `C4`, that is, edges joining `{a, b}`,
/// `{b, c}`, `{c, d}` and `{d, a}` for four distinct vertices `a, b, c, d`. Returns `c` of length
/// `m`, where `c[i]` is the number of `C4` of `G` containing the `i`-th edge.
///
/// # Complexity
/// - Time: O(n + m√m)
/// - Space: O(n + m)
///
/// # Panics
/// Panics if some `(u, v)` in `edges` satisfies `u >= n`, `v >= n` or `u = v`.
pub fn count_c4(n: usize, edges: &[(usize, usize)]) -> Vec<usize> {
    let m = edges.len();
    let mut sorted: Vec<(usize, usize, usize)> = edges
        .iter()
        .enumerate()
        .map(|(i, &(u, v))| {
            assert!(u < n, "vertex out of bounds: u={u}, n={n}");
            assert!(v < n, "vertex out of bounds: v={v}, n={n}");
            assert!(u != v, "self loop: u={u}");
            (u.min(v), u.max(v), i)
        })
        .collect();
    sorted.sort_unstable();
    let mut class = vec![0; m];
    let mut ends = vec![];
    let mut mult = vec![];
    for (k, &(u, v, i)) in sorted.iter().enumerate() {
        if k == 0 || (sorted[k - 1].0, sorted[k - 1].1) != (u, v) {
            ends.push((u, v));
            mult.push(0);
        }
        class[i] = ends.len() - 1;
        *mult.last_mut().unwrap() += 1;
    }
    let mut degree = vec![0; n];
    for &(u, v) in &ends {
        degree[u] += 1;
        degree[v] += 1;
    }
    let mut adjacency = vec![vec![]; n];
    for (e, &(u, v)) in ends.iter().enumerate() {
        adjacency[u].push((v, e));
        adjacency[v].push((u, e));
    }
    let cmp = |a: usize, b: usize| (degree[a], a) < (degree[b], b);
    let mut count = vec![0; ends.len()];
    let mut paths = vec![0; n];
    let mut touched = vec![];
    for x in 0..n {
        for &(v, e1) in &adjacency[x] {
            if cmp(v, x) {
                for &(y, e2) in &adjacency[v] {
                    if cmp(y, x) {
                        if paths[y] == 0 {
                            touched.push(y);
                        }
                        paths[y] += mult[e1] * mult[e2];
                    }
                }
            }
        }
        for &(v, e1) in &adjacency[x] {
            if cmp(v, x) {
                for &(y, e2) in &adjacency[v] {
                    if cmp(y, x) {
                        let rest = paths[y] - mult[e1] * mult[e2];
                        count[e1] += mult[e2] * rest;
                        count[e2] += mult[e1] * rest;
                    }
                }
            }
        }
        for y in touched.drain(..) {
            paths[y] = 0;
        }
    }
    class.iter().map(|&e| count[e]).collect()
}
