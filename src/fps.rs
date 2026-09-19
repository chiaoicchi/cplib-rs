use crate::algebra::multiplicative::Multiplicative;
use crate::algebra::power::pow;
use crate::algebra::{Field, Ring, RootOfUnity, Semiring};
use crate::divide_and_conquer::cdq;
use crate::periodic_function::{dft, inverse_dft};
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

impl<R: RootOfUnity<Value: PartialEq + Clone> + Default> Fps<R> {
    /// The inverse `g` of `self`, with `f g = 1`, of precision `n`.
    ///
    /// # Definition
    /// `g(0) = f(0)^{-1}` and `g(k) = -f(0)^{-1} Σ_{i=1,...,k} f(i) g(k - i)` for `k >= 1`.
    ///
    /// # Complexity
    /// - Time: O(n log n)
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `n >= 1` and `f(0) = 0`.
    pub fn inv(&self) -> Self {
        let ring = R::default();
        let n = self.precision;
        if n == 0 {
            return Self::zero(0);
        }
        assert!(
            self.value.first().is_some_and(|c| *c != ring.zero()),
            "constant term must be nonzero"
        );
        let f: Vec<R::Value> = self.coefficients().collect();
        let mut g = vec![ring.inv(&f[0])];
        let mut m = 1;
        while m < n {
            let mut e: Vec<R::Value> = f.iter().take(m << 1).cloned().collect();
            e.resize_with(m << 1, || ring.zero());
            let mut gd = g.clone();
            gd.resize_with(m << 1, || ring.zero());
            dft(&ring, &mut gd);
            dft(&ring, &mut e);
            for (x, y) in e.iter_mut().zip(&gd) {
                *x = ring.mul(x, y);
            }
            inverse_dft(&ring, &mut e);
            for x in e.iter_mut().take(m) {
                *x = ring.zero();
            }
            dft(&ring, &mut e);
            for (x, y) in e.iter_mut().zip(&gd) {
                *x = ring.mul(x, y);
            }
            inverse_dft(&ring, &mut e);
            g.extend(e[m..].iter().map(|x| ring.neg(x)));
            m <<= 1;
        }
        Self::from_vec(g, n)
    }

    /// The logarithm `g` of `self`, with `g(0) = 0` and `g' = f' / f`, of precision `n`.
    ///
    /// # Definition
    /// `log f = Σ_{k>=1} (-1)^{k+1} (f - 1)^k / k`.
    ///
    /// # Contract
    /// `1, ..., n - 1` are invertible in `R`.
    ///
    /// # Complexity
    /// - Time: O(n log n)
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `n >= 1` and `f(0) != 1`.
    pub fn log(&self) -> Self {
        let ring = R::default();
        if self.precision == 0 {
            return Self::zero(0);
        }
        assert!(
            self.value.first().is_some_and(|c| *c == ring.one()),
            "constant term must be one"
        );
        (self.derivative() * self.inv()).integral()
    }

    /// The exponential `g` of `self`, with `g(0) = 1` and `g' = f' g`, of precision `n`.
    ///
    /// # Definition
    /// `exp f = Σ_{k>=0} f^k / k!`.
    ///
    /// # Contract
    /// `1, ..., n - 1` are invertible in `R`.
    ///
    /// # Complexity
    /// - Time: O(n log n)
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `n >= 1` and `f(0) != 0`.
    pub fn exp(&self) -> Self {
        let n = self.precision;
        if n == 0 {
            return Self::zero(0);
        }
        let ring = R::default();
        assert!(
            self.value.first().is_none_or(|c| *c == ring.zero()),
            "constant term must be zero"
        );
        let f: Vec<R::Value> = self.coefficients().collect();
        let mut g = Self::one(1);
        let mut m = 1;
        while m < n {
            m = (m << 1).min(n);
            let g_lift = Self::from_vec(g.value, m);
            let f_m = Self::from_vec(f[..m].to_vec(), m);
            g = &g_lift * &(Self::one(m) - g_lift.log() + f_m);
        }
        g
    }

    /// The power `f^k`, of precision `n`.
    ///
    /// # Definition
    /// The product of `k` copies of `f` in `R[[x]]/(x^n)`, the empty product being `1`.
    ///
    /// # Contract
    /// `1`, ..., `n - 1` are invertible in `R`.
    ///
    /// # Complexity
    /// - Time: O(n log n + log k)
    /// - Space: O(n)
    pub fn pow(&self, k: u64) -> Self {
        let ring = R::default();
        let n = self.precision;
        if k == 0 {
            return Self::one(n);
        }
        let v = self.valuation();
        if v >= n || v > 0 && k >= (n as u64).div_ceil(v as u64) {
            return Self::zero(n);
        }
        let u = self.shr(v);
        let c = u.value[0].clone();
        let c_inv = ring.inv(&c);
        let u: Vec<R::Value> = u.value.iter().map(|x| ring.mul(x, &c_inv)).collect();
        let u = Self::from_vec(u, n - v);
        let k_ring = {
            let one = ring.one();
            let mut r = ring.zero();
            for i in (0..u64::BITS - k.leading_zeros()).rev() {
                r = ring.add(&r, &r);
                if k >> i & 1 == 1 {
                    r = ring.add(&r, &one);
                }
            }
            r
        };
        let l: Vec<R::Value> = u.log().value.iter().map(|x| ring.mul(x, &k_ring)).collect();
        let w = Self::from_vec(l, n - v).exp();
        let c_k = pow(&Multiplicative(R::default()), &c, k);
        let w: Vec<R::Value> = w.value.iter().map(|x| ring.mul(x, &c_k)).collect();
        let mut w = Self::from_vec(w, n - v).shl((k * v as u64) as usize);
        w.modulo(n);
        w
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

/// The solution `f` of `f(i) = step(i, Σ_{j=1,...,i} g(j) f(i - j))`, of precision `n`.
///
/// # Definition
/// `f` is the series with `f(i) = step(i, Σ_{j=1,...,i} g(j) f(i - j))` for `0 <= i < n`, the sum
/// being the coefficient of `x^i` in `f g` without the term `g(0) f(i)`. `step` is called for
/// `i = 0, 1, ..., n - 1` in this order.
///
/// # Complexity
/// - Time: O(n log^2 n), and `n` calls of `step`
/// - Space: O(n)
///
/// # Panics
/// Panics if [`poly_convolve`] panics on operands of lengths `[n/2]` and `n - 1`, that is if `R`
/// has no primitive `N`-th root of unity for the least power of two `N >= [n/2] + n - 2`.
pub fn semi_relaxed<R: RootOfUnity<Value: PartialEq + Clone> + Default>(
    g: &Fps<R>,
    mut step: impl FnMut(usize, R::Value) -> R::Value,
) -> Fps<R> {
    let ring = R::default();
    let n = g.precision;
    let g: Vec<R::Value> = g.coefficients().collect();
    let mut state: (Vec<R::Value>, Vec<R::Value>) = (
        (0..n).map(|_| ring.zero()).collect(),
        (0..n).map(|_| ring.zero()).collect(),
    );
    cdq(
        n,
        &mut state,
        |(f, acc), l, m, r| {
            let h = poly_convolve(&ring, f[l..m].to_vec(), g[1..r - l].to_vec());
            for t in m..r {
                acc[t] = ring.add(&acc[t], &h[t - l - 1]);
            }
        },
        |(f, acc), t| f[t] = step(t, acc[t].clone()),
    );
    Fps::from_vec(state.0, n)
}
