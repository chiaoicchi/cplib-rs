use crate::algebra::{Field, Ring, RootOfUnity, Semiring};
use crate::poly::Poly;
use crate::poly::poly_convolve;

/// The truncated power series ring `R[[x]]/(x^n)` over a semiring `R`.
///
/// # Definition
/// The quotient of the formal power series ring `R[[x]]` by `(x^n)`; `n` is the precision.
/// The precision of the result of an operation is the greatest `k` such that the result is
/// determined modulo `x^k` by the operands.
///
/// # Invariants
/// - `value.len() <= precision`, and `f(k) = 0` for `k >= value.len()`.
/// - The last stored coefficient is nonzero, so the zero series is the empty vector.
pub struct Fps<R: Semiring<Value: PartialEq> + Default> {
    value: Vec<R::Value>,
    precision: usize,
}

impl<R: Semiring<Value: PartialEq> + Default> Fps<R> {
    /// The image of `Σ_k value[k] x^k` in `R[[x]]/(x^precision)`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(1)
    pub fn from_vec(mut value: Vec<R::Value>, precision: usize) -> Self {
        value.truncate(precision);
        let mut fps = Self { value, precision };
        fps.normalize();
        fps
    }

    /// The image of `poly` in `R[[x]]/(x^precision)`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(1)
    pub fn from_poly(poly: Poly<R>, precision: usize) -> Self {
        Self::from_vec(poly.into_vec(), precision)
    }

    /// The zero series, the additive identity of `R[[x]]/(x^precision)`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn zero(precision: usize) -> Self {
        Self {
            value: vec![],
            precision,
        }
    }

    /// The constant series `1`, the multiplicative identity of `R[[x]]/(x^precision)`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn one(precision: usize) -> Self {
        Self::from_vec(vec![R::default().one()], precision)
    }

    /// The precision `n`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn precision(&self) -> usize {
        self.precision
    }

    /// Reduces `self` modulo `x^k`, in place, keeping the precision.
    ///
    /// # Definition
    /// `f -> Σ_{i<k} f(i)x^i`.
    ///
    /// # Complexity
    /// - Time: O(n - k)
    /// - Space: O(1)
    pub fn truncate(&mut self, k: usize) {
        self.value.truncate(k);
        self.normalize();
    }

    /// Lowers the precision of `self` to `k`, in place.
    ///
    /// # Definition
    /// The quotient map `R[[x]]/(x^n) -> R[[x]]/(x^k)`.
    ///
    /// # Complexity
    /// - Time: O(n - k)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `k > n`.
    pub fn modulo(&mut self, k: usize) {
        assert!(
            k <= self.precision,
            "k must not exceed the precision: k={k}, precision={}",
            self.precision
        );
        self.precision = k;
        self.truncate(k);
    }

    /// The formal derivative of `self`, of precision `n - 1`.
    ///
    /// # Definition
    /// The `R`-linear map `x^k -> k x^{k-1}`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `n = 0`.
    pub fn derivative(&self) -> Self {
        let ring = R::default();
        let one = ring.one();
        let mut k = ring.zero();
        Self::from_vec(
            self.value
                .iter()
                .skip(1)
                .map(|v| {
                    k = ring.add(&k, &one);
                    ring.mul(&k, v)
                })
                .collect(),
            self.precision - 1,
        )
    }

    /// The polynomial `Σ_{k<n} f(k) x^k`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn into_poly(self) -> Poly<R> {
        Poly::from_vec(self.value)
    }

    /// Restores the invariant by dropping the trailing zeros.
    ///
    /// # Complexity
    /// - Time: O(d), where `d` is the number of coefficients dropped
    /// - Space: O(1)
    fn normalize(&mut self) {
        let zero = R::default().zero();
        let len = self
            .value
            .iter()
            .rposition(|v| *v != zero)
            .map_or(0, |i| i + 1);
        self.value.truncate(len);
    }

    /// The valuation of `self`, the least `k` with `f(k) != 0`, or `n` for the zero series.
    ///
    /// # Complexity
    /// - Time: O(v)
    /// - Space: O(1)
    fn valuation(&self) -> usize {
        let zero = R::default().zero();
        self.value
            .iter()
            .position(|v| *v != zero)
            .unwrap_or(self.precision)
    }
}

impl<R: Semiring<Value: PartialEq + Clone> + Default> Fps<R> {
    /// The product `x^j f`, of precision `n + j`.
    ///
    /// # Complexity
    /// - Time: O(n + j)
    /// - Space: O(n + j)
    pub fn shl(&self, j: usize) -> Self {
        if self.value.is_empty() {
            return Self::zero(self.precision + j);
        }
        let ring = R::default();
        let mut value = Vec::with_capacity(self.value.len() + j);
        value.resize_with(j, || ring.zero());
        value.extend(self.value.iter().cloned());
        Self {
            value,
            precision: self.precision + j,
        }
    }

    /// The quotient `f / x^j`, of precision `n - j`.
    ///
    /// # Definition
    /// `f -> Σ_{k>=j} f(k) x^{k-j}.
    ///
    /// # Contract
    /// `x^j` divides `f`, that is `f(0) = ... = f(j - 1) = 0`.
    ///
    /// # Complexity
    /// - Time: O(n - j)
    /// - Space: O(n - j)
    pub fn shr(&self, j: usize) -> Self {
        Self {
            value: self
                .value
                .get(j..)
                .map_or_else(Vec::new, <[R::Value]>::to_vec),
            precision: self.precision.saturating_sub(j),
        }
    }

    /// The coefficients `f(0), ..., f(n - 1)`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn coefficients(&self) -> impl Iterator<Item = R::Value> + '_ {
        let ring = R::default();
        self.value
            .iter()
            .cloned()
            .chain(std::iter::repeat_with(move || ring.zero()))
            .take(self.precision)
    }
}

impl<R: Field<Value: PartialEq> + Default> Fps<R> {
    /// The integral of `self` vanishing at `0`, of precision `n + 1`.
    ///
    /// # Definition
    /// The `R`-linear map `x^k -> x^{k+1} / (k + 1)`.
    ///
    /// # Contract
    /// `1, ..., n` are invertible in `R`, that is the characteristic is `0` or greater than `n`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    pub fn integral(&self) -> Self {
        let ring = R::default();
        let n = self.value.len();
        if n == 0 {
            return Self::zero(self.precision + 1);
        }
        let one = ring.one();
        let minus_one = ring.neg(&one);
        let mut fact = Vec::with_capacity(n + 1);
        fact.push(ring.one());
        let mut k = ring.zero();
        for i in 1..=n {
            k = ring.add(&k, &one);
            fact.push(ring.mul(&fact[i - 1], &k));
        }
        let mut acc = ring.inv(&fact[n]);
        let mut value: Vec<R::Value> = (0..=n).map(|_| ring.zero()).collect();
        for i in (0..n).rev() {
            value[i + 1] = ring.mul(&self.value[i], &ring.mul(&fact[i], &acc));
            acc = ring.mul(&acc, &k);
            k = ring.add(&k, &minus_one);
        }
        Self {
            value,
            precision: self.precision + 1,
        }
    }
}

impl<R: Semiring<Value: PartialEq> + Default> PartialEq for Fps<R> {
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(1)
    fn eq(&self, other: &Self) -> bool {
        self.precision == other.precision && self.value == other.value
    }
}

impl<R: Semiring<Value: PartialEq + Clone> + Default> Clone for Fps<R> {
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    fn clone(&self) -> Self {
        Self {
            value: self.value.clone(),
            precision: self.precision,
        }
    }
}

impl<R: Ring<Value: PartialEq> + Default> std::ops::Neg for Fps<R> {
    type Output = Self;
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(1)
    fn neg(mut self) -> Self::Output {
        let ring = R::default();
        for v in self.value.iter_mut() {
            *v = ring.neg(v);
        }
        self
    }
}
impl<R: Ring<Value: PartialEq> + Default> std::ops::Neg for &Fps<R> {
    type Output = Fps<R>;
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    fn neg(self) -> Fps<R> {
        let ring = R::default();
        Fps {
            value: self.value.iter().map(|x| ring.neg(x)).collect(),
            precision: self.precision,
        }
    }
}

impl<R: Semiring<Value: PartialEq> + Default> std::ops::AddAssign<&Self> for Fps<R> {
    /// # Complexity
    /// - Time: O(n + m)
    /// - Space: O(1)
    fn add_assign(&mut self, rhs: &Self) {
        let ring = R::default();
        self.precision = self.precision.min(rhs.precision);
        self.value.truncate(self.precision);
        let len = rhs.value.len().min(self.precision);
        if self.value.len() < len {
            self.value.resize_with(len, || ring.zero());
        }
        for (l, r) in self.value.iter_mut().zip(rhs.value.iter()) {
            *l = ring.add(l, r);
        }
        self.normalize();
    }
}
impl<R: Ring<Value: PartialEq> + Default> std::ops::SubAssign<&Self> for Fps<R> {
    /// # Complexity
    /// - Time: O(n + m)
    /// - Space: O(1)
    fn sub_assign(&mut self, rhs: &Self) {
        let ring = R::default();
        self.precision = self.precision.min(rhs.precision);
        self.value.truncate(self.precision);
        let len = rhs.value.len().min(self.precision);
        if self.value.len() < len {
            self.value.resize_with(len, || ring.zero());
        }
        for (l, r) in self.value.iter_mut().zip(rhs.value.iter()) {
            *l = ring.add(l, &ring.neg(r));
        }
        self.normalize();
    }
}
macro_rules! forward_op_assign {
    ($($trait:ident, $method:ident, $scalars:ident);* $(;)?) => {
        $(
            impl<R: $scalars<Value: PartialEq> + Default> std::ops::$trait for Fps<R> {
                /// # Complexity
                /// - Time: O(n + m)
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

impl<R: Ring<Value: PartialEq> + RootOfUnity + Default> std::ops::MulAssign<Self> for Fps<R> {
    /// The product of precision `min(n + v(g), m + v(f))` for the valuations `v`.
    ///
    /// # Complexity
    /// - Time: O((n + m) log (n + m))
    /// - Space: O(n + m)
    fn mul_assign(&mut self, mut rhs: Self) {
        let precision = (self.precision + rhs.valuation()).min(rhs.precision + self.valuation());
        self.value.truncate(precision);
        rhs.value.truncate(precision);
        let lhs = std::mem::take(&mut self.value);
        self.value = poly_convolve(&R::default(), lhs, rhs.value);
        self.precision = precision;
        self.truncate(precision);
    }
}
impl<R: Ring<Value: PartialEq + Clone> + RootOfUnity + Default> std::ops::MulAssign<&Self>
    for Fps<R>
{
    /// # Complexity
    /// - Time: O((n + m) log (n + m))
    /// - Space: O(n + m)
    fn mul_assign(&mut self, rhs: &Self) {
        *self *= rhs.clone();
    }
}

macro_rules! forward_binop {
    ($($trait:ident, $method:ident, $assign:ident, $assign_method:ident, $scalars:ident);* $(;)?) => {
        $(
            impl<R: $scalars<Value: PartialEq> + Default> std::ops::$trait for Fps<R> {
                type Output = Self;
                fn $method(mut self, rhs: Self) -> Self {
                    std::ops::$assign::$assign_method(&mut self, rhs);
                    self
                }
            }
            impl<R: $scalars<Value: PartialEq + Clone> + Default> std::ops::$trait<&Self> for Fps<R> {
                type Output = Self;
                fn $method(mut self, rhs: &Self) -> Self {
                    std::ops::$assign::$assign_method(&mut self, rhs);
                    self
                }
            }
            impl<R: $scalars<Value: PartialEq + Clone> + Default> std::ops::$trait<Fps<R>> for &Fps<R> {
                type Output = Fps<R>;
                fn $method(self, rhs: Fps<R>) -> Fps<R> {
                    std::ops::$trait::$method(self.clone(), rhs)
                }
            }
            impl<R: $scalars<Value: PartialEq + Clone> + Default> std::ops::$trait<Self> for &Fps<R> {
                type Output = Fps<R>;
                fn $method(self, rhs: Self) -> Fps<R> {
                    std::ops::$trait::$method(self.clone(), rhs)
                }
            }
        )*
    };
}
forward_binop! {
    Add, add, AddAssign, add_assign, Semiring;
    Sub, sub, SubAssign, sub_assign, Ring;
    Mul, mul, MulAssign, mul_assign, RootOfUnity;
}
