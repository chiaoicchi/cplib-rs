use crate::algebra::RootOfUnity;
use crate::fps::elementary::fps_inv;
use crate::periodic_function::{dif, dit, twiddles};
use crate::poly::poly_convolve;

/// The subproduct tree of the points `x(0), ..., x(n - 1)`.
///
/// # Definition
/// The complete binary tree with `size` leaves, `size` the least power of two at least `n`. Leaf
/// `j` carries `t - x(j)`, with `x(j) = 0` for `j >= n`, and an inner node carries the product of
/// its two children. The root carries `t^{size - n} Π_j (t - x(j))`.
///
/// # Invariants
/// - Node `i` in heap order, `1 <= i < 2 size`, carries a monic polynomial `g_i` of degree
///   `d_i = size / 2^depth(i)`; leaf `j` is node `size + j`.
/// - `root` is the coefficients of `g_1`, of length `size + 1`.
/// - `nodes` is, for each depth from `1` to `log size` and each node `i` of that depth in order,
///   the values of `g_i` at the `2 d_i`-th roots of unity, in the order [`dif`] produces.
/// - `table` and `inv_table` are [`twiddles`] for the forward and the inverse transform.
pub struct SubproductTree<R: RootOfUnity> {
    ring: R,
    n: usize,
    size: usize,
    root: Vec<R::Value>,
    nodes: Vec<R::Value>,
    table: Vec<R::Value>,
    inv_table: Vec<R::Value>,
}

impl<R: RootOfUnity> SubproductTree<R> {
    /// Constructs the subproduct tree of `xs`.
    ///
    /// # Complexity
    /// - Time: O(n log^2 n)
    /// - Space: O(n log n)
    ///
    /// # Panics
    /// Panics if `R` has no primitive root of unity of the order the transforms need.
    pub fn new(ring: R, xs: &[R::Value]) -> Self {
        let n = xs.len();
        let size = n.next_power_of_two();
        let depth = size.trailing_zeros() as usize;
        let omega = ring
            .root_of_unity(size)
            .unwrap_or_else(|| panic!("no primitive {size}-th root of unity"));
        let inv_omega = ring.inv(&omega);
        let table = twiddles(&ring, omega, size);
        let inv_table = twiddles(&ring, inv_omega, size);

        let one = ring.one();
        let minus_one = ring.neg(&one);
        let mut nodes: Vec<R::Value> = (0..size * depth << 1).map(|_| ring.zero()).collect();
        if depth >= 1 {
            let leaves = &mut nodes[size * (depth - 1) << 1..];
            for (p, leaf) in leaves.chunks_exact_mut(2).enumerate() {
                let neg_x = if p < n { ring.neg(&xs[p]) } else { ring.zero() };
                leaf[0] = ring.add(&one, &neg_x);
                leaf[1] = ring.add(&minus_one, &neg_x);
            }
        }

        let half = ring.inv(&ring.add(&one, &one));
        let mut inv_d = ring.one();
        for l in (1..depth).rev() {
            let d = size >> l;
            inv_d = ring.mul(&inv_d, &half);
            let (upper, lower) = nodes.split_at_mut(size * l << 1);
            let parents = &mut upper[size * (l - 1) << 1..];
            let children = &lower[..size << 1];
            for (node, pair) in parents
                .chunks_exact_mut(d << 1)
                .zip(children.chunks_exact(d << 1))
            {
                let (left, right) = pair.split_at(d);
                let (lo, hi) = node.split_at_mut(d);
                for j in 0..d {
                    lo[j] = ring.mul(&left[j], &right[j]);
                    hi[j] = ring.add(&lo[j], &minus_one);
                }
                dit(&ring, hi, &inv_table);
                for (h, zeta) in hi.iter_mut().zip(&table[d..d << 1]) {
                    *h = ring.mul(&ring.mul(h, &inv_d), zeta);
                }
                dif(&ring, hi, &table);
                for h in hi.iter_mut() {
                    *h = ring.add(h, &minus_one);
                }
            }
        }

        let root = if depth == 0 {
            let neg_x = if n == 1 {
                ring.neg(&xs[0])
            } else {
                ring.zero()
            };
            vec![neg_x, ring.one()]
        } else {
            let (left, right) = nodes[..size << 1].split_at(size);
            let mut h: Vec<R::Value> = left
                .iter()
                .zip(right)
                .map(|(a, b)| ring.add(&ring.mul(a, b), &minus_one))
                .collect();
            dit(&ring, &mut h, &inv_table);
            let inv_size = ring.mul(&inv_d, &half);
            for h in h.iter_mut() {
                *h = ring.mul(h, &inv_size);
            }
            h.push(ring.one());
            h
        };

        Self {
            ring,
            n,
            size,
            root,
            nodes,
            table,
            inv_table,
        }
    }

    /// The product `Π_j (t - x(j))`, as its `n + 1` coefficients.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn product(&self) -> &[R::Value] {
        &self.root[self.size - self.n..]
    }

    /// The values `f(x(0)), ..., f(x(n - 1))` of the polynomial `Σ_k f(k) t^k`.
    ///
    /// # Complexity
    /// - Time: O(m log m + n log^2 n), where `m = f.len()`.
    /// - Space: O(m + n)
    ///
    /// # Panics
    /// Panics if `R` has no primitive root of unity of order the least power of two at least `2m`.
    pub fn evaluate(&self, f: &[R::Value]) -> Vec<R::Value> {
        let m = f.len();
        if self.n == 0 {
            return Vec::new();
        }
        if m == 0 {
            return (0..self.n).map(|_| self.ring.zero()).collect();
        }
        let zero = self.ring.zero();

        let den: Vec<R::Value> = self
            .root
            .iter()
            .rev()
            .map(|c| self.ring.add(c, &zero))
            .collect();
        let mut inv_den = fps_inv(&self.ring, &den, m);
        inv_den.reverse();
        let numer: Vec<R::Value> = f.iter().map(|c| self.ring.add(c, &zero)).collect();
        let mut cur: Vec<R::Value> = poly_convolve(&self.ring, numer, inv_den)
            .into_iter()
            .skip(m - 1)
            .take(self.size)
            .collect();
        cur.resize_with(self.size, || self.ring.zero());

        let one = self.ring.one();
        let half = self.ring.inv(&self.ring.add(&one, &one));
        let mut scale = self.ring.one();
        let mut d = self.size;
        while d > 1 {
            scale = self.ring.mul(&scale, &half);
            d >>= 1;
        }
        let mut d = self.size;
        let mut offset = 0;
        let mut next: Vec<R::Value> = (0..self.size).map(|_| self.ring.zero()).collect();
        let mut tmp: Vec<R::Value> = (0..self.size).map(|_| self.ring.zero()).collect();
        while d > 1 {
            let h = d >> 1;
            for (p, block) in cur.chunks_exact_mut(d).enumerate() {
                if p * d >= self.n {
                    break;
                }
                dif(&self.ring, block, &self.table);
                let (left, right) = self.nodes[offset + (p * d << 1)..][..d << 1].split_at(d);
                for (sibling, child) in [(right, p << 1), (left, (p << 1) + 1)] {
                    if child * h >= self.n {
                        continue;
                    }
                    let t = &mut tmp[..d];
                    for j in 0..d {
                        t[j] = self.ring.mul(&block[j], &sibling[j]);
                    }
                    dit(&self.ring, t, &self.inv_table);
                    for (o, v) in next[child * h..][..h].iter_mut().zip(&t[h..]) {
                        *o = self.ring.mul(v, &scale);
                    }
                }
            }
            std::mem::swap(&mut cur, &mut next);
            scale = self.ring.add(&scale, &scale);
            offset += self.size << 1;
            d = h;
        }
        cur.truncate(self.n);
        cur
    }

    /// The coefficients of the polynomial `f` of degree less than `n` with `f(x(j)) = y(j)`.
    ///
    /// # Definition
    /// `f = Σ_j y(j) Π_{i!=j} (t - x(i)) / Π_{i!=j} (x(j) - x(i))`, of length `n`. It is the
    /// inverse of [`evaluate`](Self::evaluate) on the polynomials of degree less than `n`.
    ///
    /// # Contract
    /// The points `x(j)` are distinct.
    ///
    /// # Complexity
    /// - Time: O(n log^2 n)
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `ys.len() != n`.
    /// Panics if `R` has no primitive root of unity of order the least power of two at least `2n`.
    pub fn interpolate(&self, ys: &[R::Value]) -> Vec<R::Value> {
        assert!(
            ys.len() == self.n,
            "length mismatch: ys.len()={}, n={}",
            ys.len(),
            self.n
        );
        if self.n == 0 {
            return Vec::new();
        }
        let zero = self.ring.zero();
        let one = self.ring.one();

        let mut k = self.ring.zero();
        let derivative: Vec<R::Value> = self.product()[1..]
            .iter()
            .map(|c| {
                k = self.ring.add(&k, &one);
                self.ring.mul(&k, c)
            })
            .collect();
        let den = self.evaluate(&derivative);

        let mut q = Vec::with_capacity(self.n);
        q.push(self.ring.mul(&den[0], &one));
        for j in 1..self.n {
            q.push(self.ring.mul(&q[j - 1], &den[j]));
        }
        let mut inv = self.ring.inv(&q[self.n - 1]);
        let mut cur: Vec<R::Value> = (0..self.size).map(|_| self.ring.zero()).collect();
        for j in (1..self.n).rev() {
            cur[j] = self.ring.mul(&ys[j], &self.ring.mul(&inv, &q[j - 1]));
            inv = self.ring.mul(&inv, &den[j]);
        }
        cur[0] = self.ring.mul(&ys[0], &inv);

        let half = self.ring.inv(&self.ring.add(&one, &one));
        let mut inv_h = self.ring.one();
        let mut h = 1;
        let mut offset = self.size * (self.size.trailing_zeros() as usize).saturating_sub(1) << 1;
        let mut next: Vec<R::Value> = (0..self.size).map(|_| self.ring.zero()).collect();
        let mut tmp: Vec<R::Value> = (0..self.size).map(|_| self.ring.zero()).collect();
        while h < self.size {
            let d = h << 1;
            for ((parent, pair), g) in next
                .chunks_exact_mut(d)
                .zip(cur.chunks_exact(d))
                .zip(self.nodes[offset..][..self.size << 1].chunks_exact(d << 1))
            {
                let (g_left, g_right) = g.split_at(d);
                for (child, sibling, first) in
                    [(&pair[..h], g_right, true), (&pair[h..], g_left, false)]
                {
                    let (lo, hi) = tmp[..d].split_at_mut(h);
                    for j in 0..h {
                        lo[j] = self.ring.add(&child[j], &zero);
                        hi[j] = self.ring.add(&child[j], &zero);
                    }
                    dit(&self.ring, hi, &self.inv_table);
                    for (v, zeta) in hi.iter_mut().zip(&self.table[h..d]) {
                        *v = self.ring.mul(&self.ring.mul(v, &inv_h), zeta);
                    }
                    dif(&self.ring, hi, &self.table);
                    for j in 0..d {
                        let term = self.ring.mul(&tmp[j], &sibling[j]);
                        parent[j] = if first {
                            term
                        } else {
                            self.ring.add(&parent[j], &term)
                        };
                    }
                }
            }
            std::mem::swap(&mut cur, &mut next);
            inv_h = self.ring.mul(&inv_h, &half);
            offset = offset.saturating_sub(self.size << 1);
            h = d;
        }

        dit(&self.ring, &mut cur, &self.inv_table);
        cur.drain(..self.size - self.n);
        for c in cur.iter_mut() {
            *c = self.ring.mul(c, &inv_h);
        }
        cur
    }
}
