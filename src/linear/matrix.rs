use crate::algebra::{One, Ring, Semiring, Zero};

/// A matrix over a semiring `R`, whose shape is determined at compile time.
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
/// The operations are taken from `R::default()`, so `R` must not carry state that distinguishes its
/// instances; a semiring whose operations are determined at run time belongs to
/// [`DynMatrix`](crate::linear::dyn_matrix::DynMatrix) instead.
///
/// # Complexity
/// - Space: O(NM)
pub struct Matrix<R: Semiring, const N: usize, const M: usize> {
    value: [[R::Value; M]; N],
    _semiring: std::marker::PhantomData<R>,
}

impl<R: Semiring + Default, const N: usize, const M: usize> Matrix<R, N, M> {
    /// Creates the zero matrix `O`, whose entries the zero element of `R`.
    ///
    /// # Complexity
    /// - Time: O(NM)
    /// - Space: O(NM)
    pub fn zero() -> Self {
        let semiring = R::default();
        Self {
            value: std::array::from_fn(|_| std::array::from_fn(|_| semiring.zero())),
            _semiring: std::marker::PhantomData,
        }
    }
}

impl<R: Semiring, const N: usize, const M: usize> Matrix<R, N, M> {
    /// Creates a matrix whose rows are `v`.
    ///
    /// # Complexity
    /// - Time: O(NM)
    /// - Space: O(NM)
    pub fn from_array(v: [[R::Value; M]; N]) -> Self {
        Self {
            value: v,
            _semiring: std::marker::PhantomData,
        }
    }

    /// Returns `true` if the matrix is square.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub const fn is_square(&self) -> bool {
        N == M
    }

    /// Returns the shape of the matrix.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub const fn shape(&self) -> (usize, usize) {
        (N, M)
    }

    /// Returns an iterator over the rows.
    ///
    /// # Complexity
    /// - Time: O(NM)
    /// - Space: O(1)
    pub fn iter(&self) -> impl Iterator<Item = &[R::Value; M]> {
        self.value.iter()
    }

    /// Returns an iterator over the rows, mutably.
    ///
    /// # Complexity
    /// - Time: O(NM)
    /// - Space: O(1)
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut [R::Value; M]> {
        self.value.iter_mut()
    }
}

impl<R: Semiring + Default, const N: usize> Matrix<R, N, N> {
    /// Creates the identity matrix `E`, whose `(i, j)` entry is `δ_ij`.
    ///
    /// # Complexity
    /// - Time: O(N^2)
    /// - Space: O(N^2)
    pub fn one() -> Self {
        let semiring = R::default();
        let value = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                if i == j {
                    semiring.one()
                } else {
                    semiring.zero()
                }
            })
        });
        Self {
            value,
            _semiring: std::marker::PhantomData,
        }
    }

    /// Raises `self` to the power of `exp`.
    ///
    /// # Complexity
    /// - Time: O(N^3 log exp)
    /// - Space: O(N^2)
    pub fn pow(&self, exp: u64) -> Self {
        let mut x = Self::one();
        for i in (0..u64::BITS - exp.leading_zeros()).rev() {
            x = &x * &x;
            if exp >> i & 1 == 1 {
                x = &x * self;
            }
        }
        x
    }
}

impl<R: Semiring, const N: usize, const M: usize> Clone for Matrix<R, N, M>
where
    R::Value: Clone,
{
    /// # Complexity
    /// - Time: O(NM)
    /// - Space: O(NM)
    fn clone(&self) -> Self {
        Self {
            value: self.value.clone(),
            _semiring: std::marker::PhantomData,
        }
    }
}

impl<R: Ring + Default, const N: usize, const M: usize> std::ops::Neg for Matrix<R, N, M> {
    type Output = Self;
    /// # Complexity
    /// - Time: O(NM)
    /// - Space: O(1)
    fn neg(mut self) -> Self::Output {
        let ring = R::default();
        for v in self.iter_mut() {
            for vi in v.iter_mut() {
                *vi = ring.neg(vi);
            }
        }
        self
    }
}
impl<R: Ring + Default, const N: usize, const M: usize> std::ops::Neg for &Matrix<R, N, M> {
    type Output = Matrix<R, N, M>;
    /// # Complexity
    /// - Time: O(NM)
    /// - Space: O(NM)
    fn neg(self) -> Matrix<R, N, M> {
        let ring = R::default();
        let value = std::array::from_fn(|i| std::array::from_fn(|j| ring.neg(&self[i][j])));
        Matrix {
            value,
            _semiring: std::marker::PhantomData,
        }
    }
}

impl<R: Semiring + Default, const N: usize, const M: usize> std::ops::Add<&Matrix<R, N, M>>
    for Matrix<R, N, M>
{
    type Output = Self;
    /// # Complexity
    /// - Time: O(NM)
    /// - Space: O(1)
    fn add(mut self, rhs: &Self) -> Self {
        let semiring = R::default();
        for (l, r) in self.iter_mut().zip(rhs.iter()) {
            for (li, ri) in l.iter_mut().zip(r.iter()) {
                *li = semiring.add(li, ri);
            }
        }
        self
    }
}
impl<R: Ring + Default, const N: usize, const M: usize> std::ops::Sub<&Matrix<R, N, M>>
    for Matrix<R, N, M>
{
    type Output = Self;
    /// # Complexity
    /// - Time: O(NM)
    /// - Space: O(1)
    fn sub(mut self, rhs: &Self) -> Self {
        let ring = R::default();
        for (l, r) in self.iter_mut().zip(rhs.iter()) {
            for (li, ri) in l.iter_mut().zip(r.iter()) {
                *li = ring.add(li, &ring.neg(ri));
            }
        }
        self
    }
}
macro_rules! forward_binop {
    ($($trait:ident, $method:ident, $scalars:ident);* $(;)?) => {
        $(
            impl<R: $scalars + Default, const N: usize, const M: usize> std::ops::$trait for Matrix<R,N, M> {
                type Output = Self;
                /// # Complexity
                /// - Time: O(NM)
                /// - Space: O(1)
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

impl<R: Semiring + Default, const N: usize, const M: usize> std::ops::Add<&Matrix<R, N, M>>
    for &Matrix<R, N, M>
{
    type Output = Matrix<R, N, M>;
    /// # Complexity
    /// - Time: O(NM)
    /// - Space: O(NM)
    fn add(self, rhs: &Matrix<R, N, M>) -> Matrix<R, N, M> {
        let semiring = R::default();
        let value =
            std::array::from_fn(|i| std::array::from_fn(|j| semiring.add(&self[i][j], &rhs[i][j])));
        Matrix {
            value,
            _semiring: std::marker::PhantomData,
        }
    }
}
impl<R: Ring + Default, const N: usize, const M: usize> std::ops::Sub<&Matrix<R, N, M>>
    for &Matrix<R, N, M>
{
    type Output = Matrix<R, N, M>;
    /// # Complexity
    /// - Time: O(NM)
    /// - Space: O(NM)
    fn sub(self, rhs: &Matrix<R, N, M>) -> Matrix<R, N, M> {
        let ring = R::default();
        let value = std::array::from_fn(|i| {
            std::array::from_fn(|j| ring.add(&self[i][j], &ring.neg(&rhs[i][j])))
        });
        Matrix {
            value,
            _semiring: std::marker::PhantomData,
        }
    }
}
macro_rules! forward_ref_binop {
    ($($trait:ident, $method:ident, $scalars:ident);* $(;)?) => {
        $(
            impl<R: $scalars + Default, const N: usize, const M: usize> std::ops::$trait<Matrix<R, N, M>> for &Matrix<R, N, M> {
                type Output = Matrix<R, N, M>;
                /// # Complexity
                /// - Time: O(NM)
                /// - Space: O(NM)
                fn $method(self, rhs: Matrix<R, N, M>) -> Matrix<R, N, M> {
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

impl<R: Semiring + Default, const N: usize, const M: usize, const L: usize>
    std::ops::Mul<&Matrix<R, M, L>> for &Matrix<R, N, M>
{
    type Output = Matrix<R, N, L>;
    /// # Complexity
    /// - Time: O(NML)
    /// - Space: O(NL)
    fn mul(self, rhs: &Matrix<R, M, L>) -> Matrix<R, N, L> {
        let semiring = R::default();
        let mut x = Matrix::zero();
        for i in 0..N {
            let xi = &mut x[i];
            let ai = &self[i];
            for (aik, bk) in ai.iter().zip(rhs.iter()) {
                for (xij, bkj) in xi.iter_mut().zip(bk) {
                    *xij = semiring.add(xij, &semiring.mul(aik, bkj));
                }
            }
        }
        x
    }
}
impl<R: Semiring + Default, const N: usize, const M: usize, const L: usize>
    std::ops::Mul<Matrix<R, M, L>> for Matrix<R, N, M>
{
    type Output = Matrix<R, N, L>;
    /// # Complexity
    /// - Time: O(NML)
    /// - Space: O(NL)
    fn mul(self, rhs: Matrix<R, M, L>) -> Matrix<R, N, L> {
        std::ops::Mul::mul(&self, &rhs)
    }
}
impl<R: Semiring + Default, const N: usize, const M: usize, const L: usize>
    std::ops::Mul<&Matrix<R, M, L>> for Matrix<R, N, M>
{
    type Output = Matrix<R, N, L>;
    /// # Complexity
    /// - Time: O(NML)
    /// - Space: O(NL)
    fn mul(self, rhs: &Matrix<R, M, L>) -> Matrix<R, N, L> {
        std::ops::Mul::mul(&self, rhs)
    }
}
impl<R: Semiring + Default, const N: usize, const M: usize, const L: usize>
    std::ops::Mul<Matrix<R, M, L>> for &Matrix<R, N, M>
{
    type Output = Matrix<R, N, L>;
    /// # Complexity
    /// - Time: O(NML)
    /// - Space: O(NL)
    fn mul(self, rhs: Matrix<R, M, L>) -> Matrix<R, N, L> {
        std::ops::Mul::mul(self, &rhs)
    }
}

impl<R: Semiring + Default, const N: usize, const M: usize> std::ops::AddAssign<&Self>
    for Matrix<R, N, M>
{
    /// # Complexity
    /// - Time: O(NM)
    /// - Space: O(1)
    fn add_assign(&mut self, rhs: &Self) {
        let semiring = R::default();
        for (l, r) in self.iter_mut().zip(rhs.iter()) {
            for (li, ri) in l.iter_mut().zip(r.iter()) {
                *li = semiring.add(li, ri);
            }
        }
    }
}
impl<R: Ring + Default, const N: usize, const M: usize> std::ops::SubAssign<&Self>
    for Matrix<R, N, M>
{
    /// # Complexity
    /// - Time: O(NM)
    /// - Space: O(1)
    fn sub_assign(&mut self, rhs: &Self) {
        let ring = R::default();
        for (l, r) in self.iter_mut().zip(rhs.iter()) {
            for (li, ri) in l.iter_mut().zip(r.iter()) {
                *li = ring.add(li, &ring.neg(ri));
            }
        }
    }
}
macro_rules! forward_op_assign {
    ($($trait:ident, $method:ident, $scalars:ident);* $(;)?) => {
        $(
            impl<R: $scalars + Default, const N: usize, const M: usize> std::ops::$trait for Matrix<R, N, M> {
                /// # Complexity
                /// - Time: O(NM)
                /// - Space: O(1)
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

impl<R: Semiring + Default, const N: usize> std::ops::MulAssign<&Self> for Matrix<R, N, N> {
    /// # Complexity
    /// - Time: O(N^3)
    /// - Space: O(N^2)
    fn mul_assign(&mut self, rhs: &Self) {
        *self = std::ops::Mul::mul(&*self, rhs);
    }
}
impl<R: Semiring + Default, const N: usize> std::ops::MulAssign<Matrix<R, N, N>>
    for Matrix<R, N, N>
{
    /// # Complexity
    /// - Time: O(N^3)
    /// - Space: O(N^2)
    fn mul_assign(&mut self, rhs: Matrix<R, N, N>) {
        *self *= &rhs;
    }
}

impl<R: Semiring, const N: usize, const M: usize> std::ops::Index<usize> for Matrix<R, N, M> {
    type Output = [R::Value; M];
    /// # Panics
    /// Panics if `index` is greater than or equal to the number of rows.
    fn index(&self, index: usize) -> &Self::Output {
        assert!(
            index < N,
            "index out of bounds: index={index}, shape={:?}",
            self.shape(),
        );
        &self.value[index]
    }
}
impl<R: Semiring, const N: usize, const M: usize> std::ops::IndexMut<usize> for Matrix<R, N, M> {
    /// # Panics
    /// Panics if `index` is greater than or equal to the number of rows.
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        assert!(
            index < N,
            "index out of bounds: index={index}, shape={:?}",
            self.shape(),
        );
        &mut self.value[index]
    }
}

impl<R: Semiring, const N: usize, const M: usize> std::cmp::PartialEq for Matrix<R, N, M>
where
    R::Value: PartialEq,
{
    fn eq(&self, rhs: &Self) -> bool {
        self.value == rhs.value
    }
}
impl<R: Semiring, const N: usize, const M: usize> std::cmp::Eq for Matrix<R, N, M> where R::Value: Eq
{}

impl<R: Semiring + Default, const N: usize, const M: usize> Zero for Matrix<R, N, M> {
    fn zero() -> Self {
        Matrix::zero()
    }
}
impl<R: Semiring + Default, const N: usize> One for Matrix<R, N, N> {
    fn one() -> Self {
        Matrix::one()
    }
}
