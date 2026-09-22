use crate::algebra::{Ring, Semiring};

/// The ring `M_n(R)` of `n x n` matrics over `R`.
///
/// # Definition
/// The `n x n` matrices over `R`, stored as in [`matrix_mul`], with the entrywise addition, the
/// product of [`matrix_mul`], the zero matrix and the identity matrix `(δ)ij)`.
///
/// # Contract
/// The values have length `n^2`.
///
/// # Complexity
/// - Time: O(n^2) for `add`, `neg`, `zero` and `one`, and O(n^3) for `mul`
/// - Space: O(n^2)
pub struct Matrix<R> {
    ring: R,
    n: usize,
}
impl<R: Semiring> Semiring for Matrix<R> {
    type Value = Vec<R::Value>;
    fn zero(&self) -> Self::Value {
        (0..self.n * self.n).map(|_| self.ring.zero()).collect()
    }
    fn one(&self) -> Self::Value {
        let mut a = self.zero();
        for i in 0..self.n {
            a[i * self.n + i] = self.ring.one();
        }
        a
    }
    fn add(&self, a: &Self::Value, b: &Self::Value) -> Self::Value {
        a.iter().zip(b).map(|(x, y)| self.ring.add(x, y)).collect()
    }
    fn mul(&self, a: &Self::Value, b: &Self::Value) -> Self::Value {
        matrix_mul(&self.ring, self.n, self.n, self.n, a, b)
    }
}
impl<R: Ring> Ring for Matrix<R> {
    fn neg(&self, a: &Self::Value) -> Self::Value {
        a.iter().map(|x| self.ring.neg(x)).collect()
    }
}

/// The product of an `n x m` matrix and an `m x l` matrix over `R`.
///
/// # Definition
/// An `n x m` matrix over `R` is a family `(a_ij)` of elements of `R` indexed by `i` in `[0, n)`
/// and `j` in `[0, m)`, stored as a vector of length `nm` with `a_ij` at index `im + j`. The
/// product is `(ab)_ij = Σ_k mul(a_ik, b_kj)`, the sum taken with `add`.
///
/// # Contract
/// `a.len() = nm` and `b.len() = ml`.
///
/// # Complexity
/// - Time: O(nml)
/// - Space: O(nl)
pub fn matrix_mul<R: Semiring>(
    ring: &R,
    n: usize,
    m: usize,
    l: usize,
    a: &[R::Value],
    b: &[R::Value],
) -> Vec<R::Value> {
    let mut c: Vec<R::Value> = (0..n * l).map(|_| ring.zero()).collect();
    for i in 0..n {
        let ci = &mut c[i * l..(i + 1) * l];
        for k in 0..m {
            let aik = &a[i * m + k];
            for (cij, bkj) in ci.iter_mut().zip(&b[k * l..(k + 1) * l]) {
                *cij = ring.add(cij, &ring.mul(aik, bkj));
            }
        }
    }
    c
}
