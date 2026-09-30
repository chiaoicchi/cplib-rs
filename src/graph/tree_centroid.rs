/// The centroid of a tree, as in [`tree_centroid`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Centroid {
    /// A single vertex.
    Vertex(usize),
    /// Two adjacent vertices, in unspecified order.
    Edge(usize, usize),
}

/// The centroid of the tree given by `edges`.
///
/// # Definition
/// The tree is on `[0, n)`, and its `i`-th edge joins `u` and `v` for `edges[i] = (u, v)`. A
/// centroid is a vertex whose removal leaves only components `C` with `2|C| <= n`. The centroids
/// are a single vertex or two adjacent vertices.
///
/// # Complexity
/// - Time: O(n)
/// - Space: O(n)
///
/// # Panics
/// Panics if `n = 0` or `n >= 2^32`.
/// Panics if some `(u, v)` in `edges` satisfies `u >= n` or `v >= n`, or `edges` is not the edge
/// set of a tree on `[0, n)`.
pub fn tree_centroid(n: usize, edges: &[(usize, usize)]) -> Centroid {
    assert!(n > 0, "not a tree: n=0");
    assert!(n < 1 << 32, "n must be less than 2^32: n={n}");
    assert_eq!(
        edges.len() + 1,
        n,
        "not a tree: number of edges={}, n={n}",
        edges.len()
    );
    let mut degree = vec![0u32; n];
    let mut xor = vec![0u32; n];
    for &(u, v) in edges {
        assert!(u < n, "vertex out of bounds: u={u}, n={n}");
        assert!(v < n, "vertex out of bounds: v={v}, n={n}");
        degree[u] += 1;
        degree[v] += 1;
        xor[u] ^= v as u32;
        xor[v] ^= u as u32;
    }

    let mut size = vec![1u32; n];
    let mut centroid: Option<usize> = None;
    let mut root = 0;
    let mut peeled = 0;
    for s in 0..n {
        let mut u = s;
        while degree[u] == 1 && peeled + 1 < n {
            let p = xor[u] as usize;
            degree[u] = 0;
            degree[p] -= 1;
            xor[p] ^= u as u32;
            size[p] += size[u];
            if (size[u] as usize) << 1 >= n && centroid.is_none_or(|c| size[u] < size[c]) {
                centroid = Some(u);
            }
            peeled += 1;
            root = p;
            u = p;
        }
    }
    assert_eq!(peeled + 1, n, "not a tree: peeled={peeled}, n={n}");

    let c = centroid.unwrap_or(root);
    if (size[c] as usize) << 1 == n {
        Centroid::Edge(c, xor[c] as usize)
    } else {
        Centroid::Vertex(c)
    }
}

/// The centroid decomposition of a tree.
///
/// # Definition
/// The tree is on `[0, n)`, and its `i`-th edge joins `u` and `v` for `edges[i] = (u, v)`, both
/// given to [`CentroidDecomposition::from_edges`]. `d(u, v)` is the number of edges on the path
/// joining `u` and `v`. Starting from the whole tree, every component that arises is split by
/// removing a centroid `c` of it, as in [`tree_centroid`]; which centroid is removed is unspecified
/// when there are two. Every vertex is removed exactly once, and `C_c` is the component from which
/// `c` is removed. The centroid tree is the tree on `[0, n)` in which the parent of `c` is the
/// vertex whose removal created `C_c`, rooted at the first vertex removed.
///
/// # Invariants
/// - `vertex[start[c]..end[c]]` is `C_c`: `c` first, then the vertices of each component of `C_c`
///   minus `c` one component after another, each in breadth-first order from `c`.
/// - For `k` in `[start[c], end[c])`, `disk[k]` is `d(c, vertex[k])`, and for `k > start[c]`,
///   `vertex[start[c] + parent[k]]` is the neighbor of `vertex[k]` on the path to `c`, and
///   `edge[k]` is the index of the edge joining them. `parent[start[c]]` and `edge[start[c]]` are
///   unused.
///
/// # Complexity
/// - Space: O(n log n)
pub struct CentroidDecomposition {
    _start: Box<[usize]>,
    _end: Box<[usize]>,
    _vertex: Box<[usize]>,
    _parent: Box<[usize]>,
    _edge: Box<[usize]>,
    _dist: Box<[usize]>,
    _ancestor_start: Box<[usize]>,
    _ancestor: Box<[(usize, usize)]>,
}

impl CentroidDecomposition {
    /// The centroid decomposition of the tree given by `edges`.
    ///
    /// # Complexity
    /// - Time: O(n log n)
    /// - Space: O(n log n)
    ///
    /// # Panics
    /// Panics if `n = 0`.
    /// Panics if some `(u, v)` in `edges` satisfies `u >= n` or `v >= n`, or `edges` is not the
    /// edge set of a tree on `[0, n)`.
    pub fn from_edges(_n: usize, _edges: &[(usize, usize)]) -> Self {
        todo!();
    }
}
