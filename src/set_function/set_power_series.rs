use crate::algebra::{Ring, Semiring};
use crate::set_function::subset::{subset_mobius, subset_zeta};

/// The ring of set power series on `[n]`.
///
/// # Definition
/// `R[x_0, ..., x_{n-1}] / (x_0^2, ..., x_{n-1}^2)`, the free `R`-module on the squarefree
/// monomials, the monomial `Π_{v∈S} x_v` being indexed by `S`. Its product is the subset
/// convolution `(fg)(S) = Σ_{T⊆S} f(T)g(S\T)`, since `x_S x_T` is `x_{S∪T}` when `S∩T = ∅` and `0`
/// otherwise. It is local, with maximal ideal `{f: f(∅) = 0}` nilpotent of index `n + 1`.
///
/// # Complexity
/// - Space: O(2^n) per value
pub struct SetPowerSeries<R> {
    ring: R,
    n: usize,
}

impl<R: Semiring> SetPowerSeries<R> {
    /// Creates the ring of set power series on `[n]`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn new(ring: R, n: usize) -> Self {
        Self { ring, n }
    }
}

impl<R: Ring<Value: Clone>> SetPowerSeries<R> {
    /// The exponential of `f`.
    ///
    /// # Definition
    /// `(exp f)(S) = Σ_{partitions of S} Π_{blocks B} f(B)`. Equivalently `Σ_k f^k / k!` when `R`
    /// has characteristic zero, a finite sum since the maximal ideal is nilpotent, and a bijection
    /// from `{f: f(∅) = 0}` onto `{f: f(∅) = 1}`.
    ///
    /// # Contract
    /// `f(∅) = 0`.
    ///
    /// # Complexity
    /// - Time: O(2^n n^2)
    /// - Space: O(2^n n)
    ///
    /// # Panics
    /// Panics if `f.len()` is not `2^n`.
    pub fn exp(&self, f: &[R::Value]) -> Vec<R::Value> {
        assert_eq!(
            f.len(),
            1 << self.n,
            "f must have length 2^n: f={}, n={}",
            f.len(),
            self.n
        );
        let mut g = self.zero();
        g[0] = self.ring.one();
        for i in 0..self.n {
            let w = 1 << i;
            let high = subset_convolve(&self.ring, f[w..w << 1].to_vec(), g[..w].to_vec());
            g[w..w << 1].clone_from_slice(&high);
        }
        g
    }

    /// The logarithm of `f`, inverting [`SetPowerSeries::exp`].
    ///
    /// # Definition
    /// The inverse of `exp` from `{f: f(∅) = 1}` onto `{f: f(∅) = 0}`.
    ///
    /// # Contract
    /// `f(∅) = 1`.
    ///
    /// # Complexity
    /// - Time: O(2^n n^2)
    /// - Space: O(2^n n)
    ///
    /// # Panics
    /// Panics if `f.len()` is not `2^n`.
    pub fn log(&self, f: &[R::Value]) -> Vec<R::Value> {
        assert_eq!(
            f.len(),
            1 << self.n,
            "f must have length 2^n: f={}, n={}",
            f.len(),
            self.n
        );
        let mut g = self.zero();
        for i in 0..self.n {
            let w = 1 << i;
            let mut a = ranked_subset_zeta(&self.ring, f[w..w << 1].to_vec());
            let b = ranked_subset_zeta(&self.ring, f[..w].to_vec());
            for s in 0..w {
                for d in 0..=i {
                    for k in 0..d {
                        let x = self.ring.mul(&a[k][s], &b[d - k][s]);
                        a[d][s] = self.ring.add(&a[d][s], &self.ring.neg(&x));
                    }
                }
            }
            g[w..w << 1].clone_from_slice(&ranked_subset_mobius(&self.ring, a));
        }
        g
    }
}

impl<R: Ring<Value: Clone>> Semiring for SetPowerSeries<R> {
    type Value = Vec<R::Value>;
    /// # Complexity
    /// - Time: O(2^n)
    /// - Space: O(2^n)
    fn zero(&self) -> Vec<R::Value> {
        (0..1 << self.n).map(|_| self.ring.zero()).collect()
    }
    /// # Complexity
    /// - Time: O(2^n)
    /// - Space: O(2^n)
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
    /// Panics if `a.len()` is not `2^n`.
    /// Panics if `b.len()` is not `2^n`.
    fn add(&self, a: &Vec<R::Value>, b: &Vec<R::Value>) -> Vec<R::Value> {
        assert_eq!(
            a.len(),
            1 << self.n,
            "a must have length 2^n: a={}, n={}",
            a.len(),
            self.n
        );
        assert_eq!(
            b.len(),
            1 << self.n,
            "b must have length 2^n: b={}, n={}",
            b.len(),
            self.n
        );
        a.iter().zip(b).map(|(x, y)| self.ring.add(x, y)).collect()
    }
    /// # Complexity
    /// - Time: O(2^n n^2)
    /// - Space: O(2^n n)
    ///
    /// # Panics
    /// Panics if `a.len()` is not `2^n`.
    /// Panics if `b.len()` is not `2^n`.
    fn mul(&self, a: &Vec<R::Value>, b: &Vec<R::Value>) -> Vec<R::Value> {
        assert_eq!(
            a.len(),
            1 << self.n,
            "a must have length 2^n: a={}, n={}",
            a.len(),
            self.n
        );
        assert_eq!(
            b.len(),
            1 << self.n,
            "b must have length 2^n: b={}, n={}",
            b.len(),
            self.n
        );
        subset_convolve(&self.ring, a.to_vec(), b.to_vec())
    }
}
impl<R: Ring<Value: Clone>> Ring for SetPowerSeries<R> {
    /// # Complexity
    /// - Time: O(2^n)
    /// - Space: O(2^n)
    ///
    /// # Panics
    /// Panics if `a.len()` is not `2^n`.
    fn neg(&self, a: &Vec<R::Value>) -> Vec<R::Value> {
        assert_eq!(
            a.len(),
            1 << self.n,
            "a must have length 2^n: a={}, n={}",
            a.len(),
            self.n,
        );
        a.iter().map(|x| self.ring.neg(x)).collect()
    }
}

/// The subset convolution of `f` and `g`.
///
/// # Definition
/// `(fg)(S) = Σ_{T⊆S} f(T) g(S\T)`, the product of [`SetPowerSeries`].
///
/// # Complexity
/// - Time: O(2^n n^2)
/// - Space: O(2^n n)
///
/// # Panics
/// Panics if `f` and `g` differ in length.
/// Panics if that length is not a power of two.
pub fn subset_convolve<R: Ring>(ring: &R, f: Vec<R::Value>, g: Vec<R::Value>) -> Vec<R::Value> {
    assert_eq!(
        f.len(),
        g.len(),
        "f and g must have the same length: f={}, g={}",
        f.len(),
        g.len()
    );
    assert!(
        f.len().is_power_of_two(),
        "length must be a power of two: len={}",
        f.len()
    );
    let n = f.len().trailing_zeros() as usize;
    let (mut f, g) = (ranked_subset_zeta(ring, f), ranked_subset_zeta(ring, g));
    for s in 0..1 << n {
        for d in (0..=n).rev() {
            let mut acc = ring.mul(&f[d][s], &g[0][s]);
            for i in (0..d).rev() {
                acc = ring.add(&acc, &ring.mul(&f[i][s], &g[d - i][s]));
            }
            f[d][s] = acc;
        }
    }
    ranked_subset_mobius(ring, f)
}

/// The ranked subset zeta transform of `f`.
///
/// # Definition
/// `layers[i][S] = Σ_{T⊆S, |T|=i} f(T)`, the subset zeta transform of `f` restricted to each rank,
/// for `i` in `[0, n]`.
///
/// # Complexity
/// - Time: O(2^n n^2)
/// - Space: O(2^n n)
///
/// # Panics
/// Panics if `f.len()` is not a power of two.
pub fn ranked_subset_zeta<R: Ring>(ring: &R, f: Vec<R::Value>) -> Vec<Vec<R::Value>> {
    assert!(
        f.len().is_power_of_two(),
        "length must be a power of two: len={}",
        f.len()
    );
    let n = f.len().trailing_zeros();
    let mut layers: Vec<Vec<R::Value>> = (0..=n)
        .map(|_| (0..f.len()).map(|_| ring.zero()).collect())
        .collect();
    for (s, x) in f.into_iter().enumerate() {
        layers[s.count_ones() as usize][s] = x;
    }
    for layer in &mut layers {
        subset_zeta(ring, layer);
    }
    layers
}

/// The diagonal of the inverse ranked subset zeta transform, inverting [`ranked_subset_zeta`].
///
/// # Definition
/// `S -> Σ_{T⊆S} (-1)^{|S| - |T|} layers[|S|][T]`, which recovers `f` from `ranked_subset_zeta(f)`.
///
/// # Complexity
/// - Time: O(2^n n^2)
/// - Space: O(2^n n)
pub fn ranked_subset_mobius<R: Ring>(ring: &R, mut layers: Vec<Vec<R::Value>>) -> Vec<R::Value> {
    for layer in &mut layers {
        subset_mobius(ring, layer);
    }
    (0..layers[0].len())
        .map(|s| std::mem::replace(&mut layers[s.count_ones() as usize][s], ring.zero()))
        .collect()
}
