use crate::algebra::or::Or;
use crate::algebra::{Ring, Semiring};
use crate::convolution::{InverseTransform, Transform};

/// The ring of set power series on `n` elements, with values `Vec<R::Value>` of length `2^n`.
///
/// # Definition
/// `R[x_1, ..., x_n] / (x_1^2, ..., x_n^2)`, the free `R`-module on the subsets of `[0, n)`
/// (the monomial `Π_{v in S} x_v` is indexed by `S`) with the product extended bilinearly from
/// `x_S x_T = x_{S or T}` if `S and T = 0` and `0` otherwise, i.e. the subset convolution
/// `(f * g)(S) = Σ_{T is subset of S} f(T) g(S\T)`. Its zero is the zero vector,
/// its one is `δ_0`, and addition is pointwise.
///
/// # Complexity
/// - Space: O(2^n) per value
pub struct SetPowerSeries<R> {
    ring: R,
    n: usize,
}

impl<R: Ring> SetPowerSeries<R> {
    /// Creates a set power series structure.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn new(ring: R, n: usize) -> Self {
        Self { ring, n }
    }

    /// Returns `2^n`, the length of a value.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        1 << self.n
    }

    /// Returns `false`, since a value always has `2^n >= 1` entries.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Returns `exp f = Σ_k f^{*k} / k!`, the sum over set partitions of the products of `f`.
    ///
    /// # Definition
    /// `(exp f)(S) = Σ_{partitions of S} Π_{blocks B} f(B)`. Splitting by the block containing the
    /// largest element `i` of `S`, `(exp f)(S) = Σ_{T subset of S\{i}} f(T or {i}) (exp f)(S\{i}\T)`,
    /// which determines `exp f` on the sets containing `i` from its values on subsets of `[0, i)`.
    ///
    /// # Contract
    /// `f(0) = 0`.
    ///
    /// # Complexity
    /// - Time: O(2^n n^2)
    /// - Space: O(2^n n)
    ///
    /// # Panics
    /// Panics if `f.len()` differs from `self`.
    pub fn exp(&self, f: &[R::Value]) -> Vec<R::Value>
    where
        R::Value: Clone,
    {
        assert!(
            f.len() == self.len(),
            "length mismatch: f={}, len={}",
            f.len(),
            self.len()
        );
        let mut g = self.zero();
        g[0] = self.ring.one();
        for i in 0..self.n {
            let w = 1 << i;
            let high = subset_convolution(&self.ring, i, &f[w..w << 1], &g[..w]);
            g[w..w << 1].clone_from_slice(&high);
        }
        g
    }

    /// Returns `log F`, the inverse of `exp`.
    ///
    /// # Definition
    /// `log` is the inverse of `exp` between `{f: f(0) = 0}` and `{F: F(0) = 1}`. Splitting by the
    /// largest element `i` of `S`, the recurrence of `exp` reads
    /// `F(S) = Σ_{T subset of S\{i}} f(T or {i}) F(S\{i}\T)`, i.e. `f` on the sets containing `i`
    /// is `F` on the sets containing `i` divided by `F` on the subsets of `[0, i)` as set power
    /// series; the division needs no inverse in `R` since `F(0) = 1`.
    ///
    /// # Contract
    /// `F(0) = 1`.
    ///
    /// # Complexity
    /// - Time: O(2^n n^2)
    /// - Space: O(2^n n)
    ///
    /// # Panics
    /// Panics if `F.len()` differs from `self`.
    pub fn log(&self, f: &[R::Value]) -> Vec<R::Value>
    where
        R::Value: Clone,
    {
        assert!(
            f.len() == self.len(),
            "length mismatch: f={}, len={}",
            f.len(),
            self.len(),
        );
        let ring = &self.ring;
        let mut g = self.zero();
        for i in 0..self.n {
            let w = 1 << i;
            let mut a = ranked_zeta(ring, i, &f[w..w << 1]);
            let b = ranked_zeta(ring, i, &f[..w]);
            for s in 0..w {
                for d in 0..=i {
                    for k in 0..d {
                        let x = ring.mul(&a[k][s], &b[d - k][s]);
                        a[d][s] = ring.add(&a[d][s], &ring.neg(&x));
                    }
                }
            }
            g[w..w << 1].clone_from_slice(&ranked_mobius(ring, i, &mut a));
        }
        g
    }
}

impl<R: Ring> Semiring for SetPowerSeries<R>
where
    R::Value: Clone,
{
    type Value = Vec<R::Value>;
    fn zero(&self) -> Vec<R::Value> {
        (0..self.len()).map(|_| self.ring.zero()).collect()
    }
    fn one(&self) -> Vec<R::Value> {
        let mut e = self.zero();
        e[0] = self.ring.one();
        e
    }
    /// # Complexity
    /// - Time: O(2^n)
    /// - Space: O(2^n)
    ///
    /// # Panics
    /// Panics if the lengths of `a`, `b` and `self` differ.
    fn add(&self, a: &Vec<R::Value>, b: &Vec<R::Value>) -> Vec<R::Value> {
        assert!(
            a.len() == self.len(),
            "length mismatch: a={}, len={}",
            a.len(),
            self.len()
        );
        assert!(
            b.len() == self.len(),
            "length mismatch: b={}, len={}",
            b.len(),
            self.len()
        );
        a.iter().zip(b).map(|(x, y)| self.ring.add(x, y)).collect()
    }
    /// # Complexity
    /// - Time: O(2^n n^2)
    /// - Space: O(2^n n)
    ///
    /// # Panics
    /// Panics if the lengths of `a`, `b` and `self` differ.
    fn mul(&self, a: &Vec<R::Value>, b: &Vec<R::Value>) -> Vec<R::Value> {
        assert!(
            a.len() == self.len(),
            "length mismatch: a={}, len={}",
            a.len(),
            self.len()
        );
        assert!(
            b.len() == self.len(),
            "length mismatch: b={}, len={}",
            b.len(),
            self.len()
        );
        subset_convolution(&self.ring, self.n, a, b)
    }
}
impl<R: Ring> Ring for SetPowerSeries<R>
where
    R::Value: Clone,
{
    /// # Complexity
    /// - Time: O(2^n)
    /// - Space: O(2^n)
    ///
    /// # Panics
    /// Panics if the length of `a` differs from `self`.
    fn neg(&self, a: &Vec<R::Value>) -> Vec<R::Value> {
        assert!(
            a.len() == self.len(),
            "length mismatch: a={}, len={}",
            a.len(),
            self.len()
        );
        a.iter().map(|x| self.ring.neg(x)).collect()
    }
}

/// Returns the subset convolution `(f * g)(S) = Σ_{T subset of S} f(T)g(S\T)` of two set power
/// series on `n` elements, by the ranked zeta transform.
///
/// # Definition
/// Split `f` by rank into `f_i(S) = [|S| = i] f(S)`, and let `Z` be the subset zeta transform,
/// `(Zf)(S) = Σ_{T subset of S} f(T)`.
/// Then `Z(f_i)(S) Z(g_j)(S) = Σ_{T, U subset of S, |T|=i, |U|=j} f(T) g(U)`, so `Z^{-1}` of it
/// counts the pairs with `T or U = S`, `|T| = i`, `|U| = j`; among them
/// `|T| + |U| = |S|` iff `T and U = 0`. Hence `(f * g)(S)` is the coefficient of rank `|S|` in
/// `Z^{-1}(Z(f)(S) Z(g)(S))`, the product being taken pointwise as polynomials in the rank.
///
/// # Complexity
/// - Time: O(2^n n^2)
/// - Space: O(2^n n)
///
/// # Panics
/// Panics if the length of `f` or `g` differs from `1 << n`.
fn subset_convolution<R: Ring>(ring: &R, n: usize, f: &[R::Value], g: &[R::Value]) -> Vec<R::Value>
where
    R::Value: Clone,
{
    assert!(
        f.len() == 1 << n,
        "length mismatch: f={}, len={}",
        f.len(),
        1 << n
    );
    assert!(
        g.len() == 1 << n,
        "length mismatch: g={}, len={}",
        g.len(),
        1 << n,
    );
    let (f, g) = (ranked_zeta(ring, n, f), ranked_zeta(ring, n, g));
    let mut h: Vec<Vec<R::Value>> = (0..=n).map(|_| vec![ring.zero(); 1 << n]).collect();
    for s in 0..1 << n {
        for i in 0..=n {
            for j in 0..=n - i {
                let x = ring.mul(&f[i][s], &g[j][s]);
                h[i + j][s] = ring.add(&h[i + j][s], &x);
            }
        }
    }
    ranked_mobius(ring, n, &mut h)
}

/// Returns `Z(f_d)` for each rank `d` in `[0, n]`, where `f_d(S) = [|S| = d] f(S)` and `Z` is the
/// subset zeta transform `(Zf)(S) = Σ_{T subset of S} f(T)`.
///
/// # Complexity
/// - Time: O(2^n n^2)
/// - Space: O(2^n n)
fn ranked_zeta<R: Ring>(ring: &R, n: usize, f: &[R::Value]) -> Vec<Vec<R::Value>>
where
    R::Value: Clone,
{
    let or = Or::new();
    let mut layers: Vec<Vec<R::Value>> = (0..=n).map(|_| vec![ring.zero(); 1 << n]).collect();
    for (s, x) in f.iter().enumerate() {
        layers[s.count_ones() as usize][s] = x.clone();
    }
    for layer in &mut layers {
        or.transform(ring, layer);
    }
    layers
}

/// Returns `S -> Z^{-1}(h_{|S|})(S)`, the subset Möbius transform of each layer read on the
/// diagonal, which recovers `f` from `ranked_zeta(f)`.
///
/// # Complexity
/// - Time: O(2^n n^2)
/// - Space: O(2^n n)
fn ranked_mobius<R: Ring>(ring: &R, n: usize, layers: &mut [Vec<R::Value>]) -> Vec<R::Value>
where
    R::Value: Clone,
{
    let or = Or::new();
    for layer in &mut *layers {
        or.inverse_transform(ring, layer);
    }
    (0usize..1 << n)
        .map(|s| layers[s.count_ones() as usize][s].clone())
        .collect()
}
