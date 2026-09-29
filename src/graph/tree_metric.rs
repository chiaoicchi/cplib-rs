use crate::algebra::Zero;

/// A diameter of the tree given by `edges` and `weight`, as its length and its vertices.
///
/// # Definition
/// The tree is on `[0, n)`, and its `i`-th edge joins `u` and `v` with weight `weight(i)` for
/// `edges[i] = (u, v)`. `d(u, v)` is the total weight of the pair joining `u` and `v`, and a
/// diameter is such a path maximizing `d(u, v)`. Returns `(d(u, v), [v_0, ..., v_k])`, where
/// `v_0 = u, ..., v_k = v` are the vertices of a diameter in order along it; which one is
/// unspecified when several exist. `weight(i)` is called exactly once for each `i`.
///
/// # Complexity
/// - Time: O(n)
/// - Space: O(n)
///
/// # Panics
/// Panics if `n = 0`.
/// Panics if some `(u, v)` in `edges` satisfies `u >= n` or `v >= n`, some `weight(i)` is less than
/// `W::zero()`, or `edges` is not the edge set of a tree on `[0, n)`.
pub fn tree_diameter<W: Copy + Ord + std::ops::Add<Output = W> + Zero>(
    _n: usize,
    _edges: &[(usize, usize)],
    mut _weight: impl FnMut(usize) -> W,
) -> (W, Vec<usize>) {
    todo!();
}

/// The center of a tree, as in [`tree_center`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Center {
    /// A single vertex.
    Vertex(usize),
    /// Two adjacent vertices, in unspecified order.
    Edge(usize, usize),
}

/// The center of the tree given by `edges` and `weight`.
///
/// # Definition
/// The tree is on `[0, n)`, and its `i`-th edge joins `u` and `v` with weight `weight(i)` for
/// `edges[i] = (u, v)`. `d(u, v)` is the total weight of the path joining `u` and `v`, and
/// `ecc(v) = max_u d(v, u)` is the eccentricity of `v`. The center is the set of the vertices
/// minimizing `ecc`, which is a single vertex or two adjacent vertices. `weight(i)` is called
/// exactly once for each `i`.
///
/// # Complexity
/// - Time: O(n)
/// - Space: O(n)
///
/// # Panics
/// Panics if `n = 0`.
/// Panics if some `(u, v)` in `edges` satisfies `u >= n` or `v >= n`, some `weight(i)` is at most
/// `W::zero()`, or `edges` is not the edge set of a tree on `[0, n)`.
pub fn tree_center<W: Copy + Ord + std::ops::Add<Output = W> + Zero>(
    _n: usize,
    _edges: &[(usize, usize)],
    mut _weight: impl FnMut(usize) -> W,
) -> Center {
    todo!();
}
