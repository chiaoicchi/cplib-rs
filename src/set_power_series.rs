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
        let (n, len, ring) = (self.n, self.len(), &self.ring);
        let or = Or::<usize>::new();
        let ranked = |f: &[R::Value]| -> Vec<Vec<R::Value>> {
            let mut layers: Vec<Vec<R::Value>> = (0..=n).map(|_| self.zero()).collect();
            for (s, x) in f.iter().enumerate() {
                layers[s.count_ones() as usize][s] = x.clone();
            }
            for layer in &mut layers {
                or.transform(ring, layer);
            }
            layers
        };
        let (fa, fb) = (ranked(a), ranked(b));
        let mut h: Vec<Vec<R::Value>> = (0..=n).map(|_| self.zero()).collect();
        for s in 0..len {
            for i in 0..=n {
                for j in 0..=n - i {
                    let x = ring.mul(&fa[i][s], &fb[j][s]);
                    h[i + j][s] = ring.add(&h[i + j][s], &x);
                }
            }
        }
        for layer in &mut h {
            or.inverse_transform(ring, layer);
        }
        (0..len)
            .map(|s| h[s.count_ones() as usize][s].clone())
            .collect()
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
