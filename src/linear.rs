use crate::algebra::{Ring, Semiring};

/// A matrix over a semiring `R`, whose shape is determined at run time.
///
/// # Definition
/// An `n x m` matrix over `R` is a family `(a_ij)` of elements of `R` indexed by `0 <= i < n`
/// and `0 <= j < m`.
/// Addition is entrywise, `(a + b)_ij = add(a_ij, b_ij)`, and multiplication of an `n x m` matrix
/// by an `m x l` matrix is `(ab)_il = ∑_j mul(a_ij, b_jl)`, the sum taken with `add`.
/// The `n x n` matrices form a semiring under these operations, with the zero matrix `O` and the
/// identity matrix `E`.
///
/// # Contract
/// Both operands of a binary operation must hold the same `semiring`.
///
/// # Invariants
/// - `n > 0` and `m > 0`.
/// - `n * m = value.len()`.
///
/// # Complexity
/// - Space: O(nm)
pub struct DynMatrix<R: Semiring> {
    semiring: R,
    n: usize,
    m: usize,
    value: Vec<R::Value>,
}

impl<R: Semiring> DynMatrix<R> {
    /// Creates the zero dynanic matrix `O`, whose entries are all `semiring.zero()`.
    ///
    /// # Complexity
    /// - Time: O(nm)
    /// - Space: O(nm)
    ///
    /// # Panics
    /// Panics if `n == 0`.
    /// Panics if `m == 0`.
    pub fn zero(semiring: R, n: usize, m: usize) -> Self {
        assert!(n > 0, "n must be greater than 0");
        assert!(m > 0, "m must be greater than 0");
        Self {
            value: (0..n * m).map(|_| semiring.zero()).collect(),
            semiring,
            n,
            m,
        }
    }

    /// Creates the identity dynamic matrix `E`, whose `(i, j)` entry is `δ_ij`.
    ///
    /// # Complexity
    /// - Time: O(n^2)
    /// - Space: O(n^2)
    ///
    /// # Panics
    /// Panics if `n == 0`.
    pub fn one(semiring: R, n: usize) -> Self {
        assert!(n > 0, "n must be greater than 0");
        let mut value: Vec<R::Value> = (0..n * n).map(|_| semiring.zero()).collect();
        for i in 0..n {
            value[i * n + i] = semiring.one();
        }
        Self {
            semiring,
            n,
            m: n,
            value,
        }
    }

    /// Creates a matrix whose rows are `v`.
    ///
    /// # Complexity
    /// - Time: O(nm)
    /// - Space: O(nm)
    ///
    /// # Panics
    /// Panics if `v` is empty.
    /// Panics if any row of `v` is empty.
    /// Panics if any row of `v` is not the same length.
    pub fn from_vec(semiring: R, v: Vec<Vec<R::Value>>) -> Self {
        assert!(!v.is_empty(), "`v` must not be empty");
        assert!(!v[0].is_empty(), "row of `v` must not be empty");
        assert!(
            v.iter().all(|vi| vi.len() == v[0].len()),
            "all rows of `v` must have the same length",
        );
        Self {
            semiring,
            n: v.len(),
            m: v[0].len(),
            value: v.into_iter().flatten().collect(),
        }
    }

    /// Returns `true` if the matrix is square.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_square(&self) -> bool {
        self.n == self.m
    }

    /// Returns the shape of the matrix.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn shape(&self) -> (usize, usize) {
        (self.n, self.m)
    }

    /// Returns an iterator over the rows.
    ///
    /// # Complexity
    /// - Time: O(nm)
    /// - Space: O(1)
    pub fn iter(&self) -> impl Iterator<Item = &[R::Value]> {
        self.value.chunks_exact(self.m)
    }

    /// Returns an iterator over the rows, mutably.
    ///
    /// # Complexity
    /// - Time: O(nm)
    /// - Space: O(1)
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut [R::Value]> {
        self.value.chunks_exact_mut(self.m)
    }
}

impl<R: Semiring + Clone> DynMatrix<R> {
    /// Raises `self` to the power of `exp`.
    ///
    /// # Complexity
    /// - Time: O(n^3 log exp)
    /// - Space: O(n^2)
    ///
    /// # Panics
    /// Panics if `self` is not square.
    pub fn pow(&self, exp: u64) -> Self {
        assert!(self.is_square(), "not square: shape={:?}", self.shape());
        let mut x = Self::one(self.semiring.clone(), self.n);
        for i in (0..u64::BITS - exp.leading_zeros()).rev() {
            x = &x * &x;
            if exp >> i & 1 == 1 {
                x = &x * self;
            }
        }
        x
    }
}

impl<R: Semiring + Clone> Clone for DynMatrix<R>
where
    R::Value: Clone,
{
    /// # Complexity
    /// - Time: O(nm)
    /// - Space: O(nm)
    fn clone(&self) -> Self {
        Self {
            semiring: self.semiring.clone(),
            n: self.n,
            m: self.m,
            value: self.value.clone(),
        }
    }
}

impl<R: Ring> std::ops::Neg for DynMatrix<R> {
    type Output = Self;
    /// # Complexity
    /// - Time: O(nm)
    /// - Space: O(1)
    fn neg(mut self) -> Self::Output {
        for v in self.value.iter_mut() {
            *v = self.semiring.neg(v);
        }
        self
    }
}
impl<R: Ring + Clone> std::ops::Neg for &DynMatrix<R> {
    type Output = DynMatrix<R>;
    /// # Complexity
    /// - Time: O(nm)
    /// - Space: O(nm)
    fn neg(self) -> DynMatrix<R> {
        DynMatrix {
            semiring: self.semiring.clone(),
            n: self.n,
            m: self.m,
            value: self.value.iter().map(|x| self.semiring.neg(x)).collect(),
        }
    }
}

impl<R: Semiring> std::ops::Add<&DynMatrix<R>> for DynMatrix<R> {
    type Output = Self;
    /// # Complexity
    /// - Time: O(nm)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if the shapes of `self` and `rhs` differ.
    fn add(mut self, rhs: &Self) -> Self {
        assert!(
            self.shape() == rhs.shape(),
            "shape mismatch: lhs:{:?}, rhs:{:?}",
            self.shape(),
            rhs.shape()
        );
        for (l, r) in self.value.iter_mut().zip(rhs.value.iter()) {
            *l = self.semiring.add(l, r);
        }
        self
    }
}
impl<R: Ring> std::ops::Sub<&DynMatrix<R>> for DynMatrix<R> {
    type Output = Self;
    /// # Complexity
    /// - Time: O(nm)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if the shapes of `self` and `rhs` differ.
    fn sub(mut self, rhs: &Self) -> Self {
        assert!(
            self.shape() == rhs.shape(),
            "shape mismatch: lhs:{:?}, rhs:{:?}",
            self.shape(),
            rhs.shape(),
        );
        for (l, r) in self.value.iter_mut().zip(rhs.value.iter()) {
            *l = self.semiring.add(l, &self.semiring.neg(r));
        }
        self
    }
}
macro_rules! forward_binop {
    ($($trait:ident, $method:ident, $scalars:ident);* $(;)?) => {
        $(
            impl<R: $scalars> std::ops::$trait for DynMatrix<R> {
                type Output = Self;
                /// # Complexity
                /// - Time: O(nm)
                /// - Space: O(1)
                ///
                /// # Panics
                /// Panics if the shapes of `self` and `rhs` differ.
                fn $method(self, rhs: Self) -> Self {
                    std::ops::$trait::$method(self, &rhs)
                }
            }
        )*
    };
}
forward_binop! {
    Add, add, Semiring;
    Sub, sub, Ring;
}

impl<R: Semiring + Clone> std::ops::Add<&DynMatrix<R>> for &DynMatrix<R> {
    type Output = DynMatrix<R>;
    /// # Complexity
    /// - Time: O(nm)
    /// - Space: O(nm)
    ///
    /// # Panics
    /// Panics if the shapes of `self` and `rhs` differ.
    fn add(self, rhs: &DynMatrix<R>) -> DynMatrix<R> {
        assert!(
            self.shape() == rhs.shape(),
            "shape mismatch: lhs:{:?}, rhs:{:?}",
            self.shape(),
            rhs.shape()
        );
        DynMatrix {
            semiring: self.semiring.clone(),
            n: self.n,
            m: self.m,
            value: self
                .value
                .iter()
                .zip(&rhs.value)
                .map(|(l, r)| self.semiring.add(l, r))
                .collect(),
        }
    }
}
impl<R: Ring + Clone> std::ops::Sub<&DynMatrix<R>> for &DynMatrix<R> {
    type Output = DynMatrix<R>;
    /// # Complexity
    /// - Time: O(nm)
    /// - Space: O(nm)
    ///
    /// # Panics
    /// Panics if the shapes of `self` and `rhs` differ.
    fn sub(self, rhs: &DynMatrix<R>) -> DynMatrix<R> {
        assert!(
            self.shape() == rhs.shape(),
            "shape mismatch: lhs:{:?}, rhs:{:?}",
            self.shape(),
            rhs.shape(),
        );
        DynMatrix {
            semiring: self.semiring.clone(),
            n: self.n,
            m: self.m,
            value: self
                .value
                .iter()
                .zip(&rhs.value)
                .map(|(l, r)| self.semiring.add(l, &self.semiring.neg(r)))
                .collect(),
        }
    }
}
macro_rules! forward_ref_binop {
    ($($trait:ident, $method:ident, $scalars:ident);* $(;)?) => {
        $(
            impl<R: $scalars + Clone> std::ops::$trait<DynMatrix<R>> for &DynMatrix<R> {
                type Output = DynMatrix<R>;
                /// # Complexity
                /// - Time: O(nm)
                /// - Space: O(nm)
                ///
                /// # Panics
                /// Panics if the shapes of `self` and `rhs` differ.
                fn $method(self, rhs: DynMatrix<R>) -> DynMatrix<R> {
                    std::ops::$trait::$method(self, &rhs)
                }
            }
        )*
    };
}
forward_ref_binop! {
    Add, add, Semiring;
    Sub, sub, Ring;
}

impl<R: Semiring + Clone> std::ops::Mul<&DynMatrix<R>> for &DynMatrix<R> {
    type Output = DynMatrix<R>;
    /// # Complexity
    /// - Time: O(nml)
    /// - Space: O(nl)
    ///
    /// where `self.shape() = (n, m)`, `rhs.shape() = (m, l)`.
    ///
    /// # Panics
    /// Panics if the number of columns of `self` differs from the number of rows of `rhs`.
    fn mul(self, rhs: &DynMatrix<R>) -> DynMatrix<R> {
        assert!(
            self.m == rhs.n,
            "shape mismatch for multiplication: lhs:{:?}, rhs:{:?}",
            self.shape(),
            rhs.shape(),
        );
        let mut x = DynMatrix::zero(self.semiring.clone(), self.n, rhs.m);
        for i in 0..self.n {
            let xi = &mut x[i];
            let ai = &self[i];
            for (aik, bk) in ai.iter().zip(rhs.iter()) {
                for (xij, bkj) in xi.iter_mut().zip(bk) {
                    *xij = self.semiring.add(xij, &self.semiring.mul(aik, bkj));
                }
            }
        }
        x
    }
}
impl<R: Semiring + Clone> std::ops::Mul for DynMatrix<R> {
    type Output = Self;
    /// # Complexity
    /// - Time: O(nml)
    /// - Space: O(nl)
    ///
    /// where `self.shape() = (n, m)`, `rhs.shape() = (m, l)`.
    ///
    /// # Panics
    /// Panics if the number of columns of `self` differs from the number of rows of `rhs`.
    fn mul(self, rhs: Self) -> Self {
        std::ops::Mul::mul(&self, &rhs)
    }
}
impl<R: Semiring + Clone> std::ops::Mul<&DynMatrix<R>> for DynMatrix<R> {
    type Output = Self;
    /// # Complexity
    /// - Time: O(nml)
    /// - Space: O(nl)
    ///
    /// where `self.shape() = (n, m)`, `rhs.shape() = (m, l)`.
    ///
    /// # Panics
    /// Panics if the number of columns of `self` differs from the number of rows of `rhs`.
    fn mul(self, rhs: &DynMatrix<R>) -> Self {
        std::ops::Mul::mul(&self, rhs)
    }
}
impl<R: Semiring + Clone> std::ops::Mul<DynMatrix<R>> for &DynMatrix<R> {
    type Output = DynMatrix<R>;
    /// # Complexity
    /// - Time: O(nml)
    /// - Space: O(nl)
    ///
    /// where `self.shape() = (n, m)`, `rhs.shape() = (m, l)`.
    ///
    /// # Panics
    /// Panics if the number of columns of `self` differs from the number of rows of `rhs`.
    fn mul(self, rhs: DynMatrix<R>) -> DynMatrix<R> {
        std::ops::Mul::mul(self, &rhs)
    }
}

impl<R: Semiring> std::ops::AddAssign<&Self> for DynMatrix<R> {
    /// # Complexity
    /// - Time: O(nm)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if the shapes of `self` and `rhs` differ.
    fn add_assign(&mut self, rhs: &Self) {
        assert!(
            self.shape() == rhs.shape(),
            "shape mismatch: lhs:{:?}, rhs:{:?}",
            self.shape(),
            rhs.shape(),
        );
        for (l, r) in self.value.iter_mut().zip(rhs.value.iter()) {
            *l = self.semiring.add(l, r);
        }
    }
}
impl<R: Ring> std::ops::SubAssign<&Self> for DynMatrix<R> {
    /// # Complexity
    /// - Time: O(nm)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if the shapes of `self` and `rhs` differ.
    fn sub_assign(&mut self, rhs: &Self) {
        assert!(
            self.shape() == rhs.shape(),
            "shape mismatch: lhs:{:?}, rhs:{:?}",
            self.shape(),
            rhs.shape(),
        );
        for (l, r) in self.value.iter_mut().zip(rhs.value.iter()) {
            *l = self.semiring.add(l, &self.semiring.neg(r));
        }
    }
}
macro_rules! forward_op_assign {
    ($($trait:ident, $method:ident, $scalars:ident);* $(;)?) => {
        $(
            impl<R: $scalars> std::ops::$trait for DynMatrix<R> {
                /// # Complexity
                /// - Time: O(nm)
                /// - Space: O(1)
                ///
                /// # Panics
                /// Panics if the shapes of `self` and `rhs` differ.
                fn $method(&mut self, rhs: Self) {
                    std::ops::$trait::$method(self, &rhs);
                }
            }
        )*
    };
}
forward_op_assign! {
    AddAssign, add_assign, Semiring;
    SubAssign, sub_assign, Ring;
}

impl<R: Semiring + Clone> std::ops::MulAssign<&Self> for DynMatrix<R> {
    /// # Complexity
    /// - Time: O(nml)
    /// - Space: O(nl)
    ///
    /// where `self.shape() = (n, m)`, `rhs.shape() = (m, l)`.
    ///
    /// # Panics
    /// Panics if the number of columns of `self` differs from the number of rows of `rhs`.
    fn mul_assign(&mut self, rhs: &Self) {
        *self = std::ops::Mul::mul(&*self, rhs);
    }
}
impl<R: Semiring + Clone> std::ops::MulAssign for DynMatrix<R> {
    /// # Complexity
    /// - Time: O(nml)
    /// - Space: O(nl)
    ///
    /// where `self.shape() = (n, m)`, `rhs.shape() = (m, l)`.
    ///
    /// # Panics
    /// Panics if the number of columns of `self` differs from the number of rows of `rhs`.
    fn mul_assign(&mut self, rhs: Self) {
        *self *= &rhs;
    }
}

impl<R: Semiring> std::ops::Index<usize> for DynMatrix<R> {
    type Output = [R::Value];
    /// # Panics
    /// Panics if `index` is greater than or equal to the number of rows.
    fn index(&self, index: usize) -> &Self::Output {
        assert!(
            index < self.n,
            "index out of bounds: index={index}, shape={:?}",
            self.shape(),
        );
        &self.value[index * self.m..(index + 1) * self.m]
    }
}
impl<R: Semiring> std::ops::IndexMut<usize> for DynMatrix<R> {
    /// # Panics
    /// Panics if `index` is greater than or equal to the number of rows.
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        assert!(
            index < self.n,
            "index out of bounds: index={index}, shape={:?}",
            self.shape(),
        );
        &mut self.value[index * self.m..(index + 1) * self.m]
    }
}
