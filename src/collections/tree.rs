/// A rooted tree with a heavy-light DFS order.
///
/// # Definition
/// Let `T` be a tree on `[0, n)` rooted at `root`. The heavy child of a vertex is a child of
/// maximum subtree size, and an edge to a heavy child is heavy; other edges are light. The maximal
/// paths of heavy edges are the heavy paths, and the shallowest vertex of a heavy path is its head.
/// `order` is the preorder of the DFS that visits the heavy child first, and `index` is its
/// inverse.
///
/// A path from `u` to `v` is returned as a sequence of `Segment`s of `order`, in the order they are
/// traversed from `u` to `v`: `Up(r)` is traversed in decreasing order of `order` (towards the
/// root), `Down(r)` in increasing order. `edge_path` omits `lca(u, v)`, so that, if the value of
/// the edge `(parent[c], c)` is stored at position `index[c]`, the segments cover exactly the edges
/// of the path.
///
/// # Contract
/// `adjacency` is the adjacency list of an undirected tree on `[0, n)`, and `root < n`.
///
/// # Invariants
/// - `order` is the heavy-child-first preorder from `root`, and `index[order[i]] = i`.
/// - `parent[root] = root`, `depth[root] = 0`, and `head[v]` is the head of the heavy path of `v`.
/// - `out[v] = index[v] + size[v]`, so the subtree of `v` is `order[index[v]..out[v]]`.
///
/// # Complexity
/// - Space: O(n)
pub struct Tree {
    parent: Box<[usize]>,
    depth: Box<[usize]>,
    head: Box<[usize]>,
    index: Box<[usize]>,
    out: Box<[usize]>,
    order: Box<[usize]>,
}

/// A segment of a path in a Tree, as an interval of the DFS order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Segment {
    /// Traversed in decreasing order, i.e. towards the root.
    Up(std::ops::Range<usize>),
    /// Traversed in increasing order, i.e. away from the root.
    Down(std::ops::Range<usize>),
}

impl Tree {
    /// Builds the tree rooted at `root`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `root >= n`.
    /// Panics if `adjacency` is not a tree.
    pub fn from_adjacency(adjacency: &[Vec<usize>], root: usize) -> Self {
        let n = adjacency.len();
        assert!(root < n, "root out of bounds: root={root}, n={n}");
        let degree: usize = adjacency.iter().map(Vec::len).sum();
        assert!(
            degree == 2 * (n - 1),
            "not a tree: sum of degrees={degree}, n={n}"
        );
        let mut parent = vec![root; n];
        let mut depth = vec![0; n];
        let mut seen = vec![false; n];
        seen[root] = true;
        let mut preorder = Vec::with_capacity(n);
        let mut stack = vec![root];
        while let Some(u) = stack.pop() {
            preorder.push(u);
            for &v in &adjacency[u] {
                if v == parent[u] {
                    continue;
                }
                assert!(!seen[v], "not a tree: vertex {v} is reached twice");
                seen[v] = true;
                parent[v] = u;
                depth[v] = depth[u] + 1;
                stack.push(v);
            }
        }
        assert!(
            preorder.len() == n,
            "not a tree: reachable={}, n={n}",
            preorder.len()
        );

        let mut size = vec![1; n];
        let mut heavy = vec![usize::MAX; n];
        for &v in preorder.iter().rev().take(n - 1) {
            let p = parent[v];
            size[p] += size[v];
            if heavy[p] == usize::MAX || size[heavy[p]] < size[v] {
                heavy[p] = v;
            }
        }

        let mut head = vec![root; n];
        let mut index = vec![0; n];
        let mut order = Vec::with_capacity(n);
        let mut stack = vec![root];
        while let Some(u) = stack.pop() {
            index[u] = order.len();
            order.push(u);
            for &v in &adjacency[u] {
                if v != parent[u] && v != heavy[u] {
                    head[v] = v;
                    stack.push(v);
                }
            }
            if heavy[u] != usize::MAX {
                head[heavy[u]] = head[u];
                stack.push(heavy[u]);
            }
        }
        let out = (0..n).map(|v| index[v] + size[v]).collect();
        Self {
            parent: parent.into(),
            depth: depth.into(),
            head: head.into(),
            index: index.into(),
            out,
            order: order.into(),
        }
    }

    /// Returns the parent of `v`, which is `root` itself for `v = root`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `v >= n`.
    pub fn parent(&self, v: usize) -> usize {
        assert!(
            v < self.len(),
            "index out of bounds: v={v}, len={}",
            self.len()
        );
        self.parent[v]
    }

    /// Returns the depth of `v`, the number of edges from `root`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `v >= n`.
    pub fn depth(&self, v: usize) -> usize {
        assert!(
            v < self.len(),
            "index out of bounds: v={v}, len={}",
            self.len()
        );
        self.depth[v]
    }

    /// Returns the position of `v` in `order`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `v >= n`.
    pub fn index(&self, v: usize) -> usize {
        assert!(
            v < self.len(),
            "index out of bounds: v={v}, len={}",
            self.len()
        );
        self.index[v]
    }

    /// Returns the vertex at position `i` of `order`, the inverse of `index`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i >= n`.
    pub fn vertex(&self, i: usize) -> usize {
        assert!(
            i < self.len(),
            "index out of bounds: i={i}, len={}",
            self.len()
        );
        self.order[i]
    }

    /// Returns the interval of `order` occupied by the subtree of `v`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `v >= n`.
    pub fn subtree(&self, v: usize) -> std::ops::Range<usize> {
        assert!(
            v < self.len(),
            "index out of bounds: v={v}, len={}",
            self.len()
        );
        self.index[v]..self.out[v]
    }

    /// Returns `true` if `u` is an ancestor of `v`, including `u = v`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `u >= n` or `v >= n`.
    pub fn is_ancestor(&self, u: usize, v: usize) -> bool {
        assert!(
            u < self.len(),
            "index out of bounds: u={u}, len={}",
            self.len()
        );
        assert!(
            v < self.len(),
            "index out of bounds: v={v}, len={}",
            self.len()
        );
        self.subtree(u).contains(&self.index[v])
    }

    /// Returns the `k`-th ancestor of `v`, and `None` if `k > depth[v]`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `v >= n`.
    pub fn ancestor(&self, mut v: usize, mut k: usize) -> Option<usize> {
        assert!(
            v < self.len(),
            "index out of bounds: v={v}, len={}",
            self.len()
        );
        if k > self.depth[v] {
            return None;
        }
        loop {
            let h = self.head[v];
            let d = self.index[v] - self.index[h];
            if k <= d {
                return Some(self.order[self.index[v] - k]);
            }
            k -= d + 1;
            v = self.parent[h];
        }
    }

    /// Returns the lowest common ancestor of `u` and `v`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `u >= n` or `v >= n`.
    pub fn lca(&self, mut u: usize, mut v: usize) -> usize {
        assert!(
            u < self.len(),
            "index out of bounds: u={u}, len={}",
            self.len()
        );
        assert!(
            v < self.len(),
            "index out of bounds: v={v}, len={}",
            self.len()
        );
        while self.head[u] != self.head[v] {
            if self.depth[self.head[u]] > self.depth[self.head[v]] {
                u = self.parent[self.head[u]];
            } else {
                v = self.parent[self.head[v]];
            }
        }
        if self.depth[u] < self.depth[v] { u } else { v }
    }

    /// Returns the number of edges on the path from `u` to `v`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `u >= n` or `v >= n`.
    pub fn dist(&self, u: usize, v: usize) -> usize {
        assert!(
            u < self.len(),
            "index out of bounds: u={u}, len={}",
            self.len()
        );
        assert!(
            v < self.len(),
            "index out of bounds: v={v}, len={}",
            self.len()
        );
        self.depth[u] + self.depth[v] - (self.depth[self.lca(u, v)] << 1)
    }

    /// Returns the vertex `k` edges from `u` on the path to `v`, and `None` if `k > dist(u, v)`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `u >= n` or `v >= n`.
    pub fn jump(&self, u: usize, v: usize, k: usize) -> Option<usize> {
        assert!(
            u < self.len(),
            "index out of bounds: u={u}, len={}",
            self.len()
        );
        assert!(
            v < self.len(),
            "index out of bounds: v={v}, len={}",
            self.len()
        );
        let w = self.lca(u, v);
        let up = self.depth[u] - self.depth[w];
        let d = up + self.depth[v] - self.depth[w];
        if k <= up {
            self.ancestor(u, k)
        } else if k <= d {
            self.ancestor(v, d - k)
        } else {
            None
        }
    }

    /// Returns the segments covering the vertices of the path from `u` to `v`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(log n)
    ///
    /// # Panics
    /// Panics if `u >= n` or `v >= n`.
    pub fn path(&self, mut u: usize, mut v: usize) -> Vec<Segment> {
        assert!(
            u < self.len(),
            "index out of bounds: u={u}, len={}",
            self.len()
        );
        assert!(
            v < self.len(),
            "index out of bounds: v={v}, len={}",
            self.len()
        );
        let mut up = Vec::new();
        let mut down = Vec::new();
        while self.head[u] != self.head[v] {
            if self.depth[self.head[u]] >= self.depth[self.head[v]] {
                let h = self.head[u];
                up.push(Segment::Up(self.index[h]..self.index[u] + 1));
                u = self.parent[h];
            } else {
                let h = self.head[v];
                down.push(Segment::Down(self.index[h]..self.index[v] + 1));
                v = self.parent[h];
            }
        }
        let (i, j) = (self.index[u], self.index[v]);
        if i >= j {
            up.push(Segment::Up(j..i + 1));
        } else {
            down.push(Segment::Down(i..j + 1));
        }
        up.extend(down.into_iter().rev());
        up
    }

    /// Returns the segments covering the vertices of the path from `u` to `v` except `lca(u, v)`.
    ///
    /// # Complexity
    /// - Time: O(log n)
    /// - Space: O(log n)
    ///
    /// # Panics
    /// Panics if `u >= n` or `v >= n`.
    pub fn edge_path(&self, mut u: usize, mut v: usize) -> Vec<Segment> {
        assert!(
            u < self.len(),
            "index out of bounds: u={u}, len={}",
            self.len()
        );
        assert!(
            v < self.len(),
            "index out of bounds: v={v}, len={}",
            self.len()
        );
        let mut up = Vec::new();
        let mut down = Vec::new();
        while self.head[u] != self.head[v] {
            if self.depth[self.head[u]] >= self.depth[self.head[v]] {
                let h = self.head[u];
                up.push(Segment::Up(self.index[h]..self.index[u] + 1));
                u = self.parent[h];
            } else {
                let h = self.head[v];
                down.push(Segment::Down(self.index[h]..self.index[v] + 1));
                v = self.parent[h];
            }
        }
        let (i, j) = (self.index[u], self.index[v]);
        if i > j {
            up.push(Segment::Up(j + 1..i + 1));
        } else if j > i {
            down.push(Segment::Down(i + 1..j + 1));
        }
        up.extend(down.into_iter().rev());
        up
    }

    /// Returns the number of vertices.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.order.len()
    }

    /// Returns `true` if the tree has no vertices.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        self.order.is_empty()
    }
}
