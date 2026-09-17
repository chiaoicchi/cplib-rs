use crate::algebra::additive::Additive;
use crate::algebra::canonical::Canonical;
use crate::algebra::{Field, Ring, RootOfUnity, Semiring};
use crate::convolution::Convolution;

/// The polynomial ring `R[x]` over a semiring `R`.
///
/// # Definition
/// `R[x]` is the free `R`-module on the monomials `x^k` for `k >= 0`, with the product determined
/// by `x^i x^j = x^{i+j}` and extended bilinearly, so `(fg)(k) = Σ_{i+j=k} f(i)g(j)`. It is the
/// monoid algebra `R[(N, +)]`, and it is the free `R`-algebra on one generator: a homomorphism out
/// of `R[x]` is determined by the image of `x`, which is why evaluation exists.
///
/// A polynomial has finitely many nonzero coefficients, so it is stored as `f(0), ..., f(deg)`.
/// This is what separates `R[x]` from the formal power series ring `R[[x]]`, whose coefficient
/// sequences are arbitrary. The separation cuts both ways: `inv`, `log` and `exp` need the `x`-adic
/// completeness of `R[[x]]` and are not operations here, while evaluation, division with remainder
/// and `taylor_shift` need the degree and are not operations there.
///
/// # Invariants
/// - The last stored coefficient is nonzero, so the zero polynomial is the empty vector and the
///   stored length is `deg + 1`. Every operation restores this; `mul`, `shl`, `shr` and `neg`
///   preserve it without a scan, since none of them can cancel the leading coefficient.
pub struct Poly<R: Semiring<Value: PartialEq> + Default>(Vec<R::Value>);

impl<R: Semiring<Value: PartialEq> + Default> Poly<R> {
    /// Constructs `Σ value[k] x^k`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(1)
    pub fn from_vec(value: Vec<R::Value>) -> Self {
        let mut poly = Self(value);
        poly.normalize();
        poly
    }

    /// The zero polynomial.
    ///
    /// # Definition
    /// The additive identity of `R[x]`. It has no nonzero coefficient, hence no degree.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn zero() -> Self {
        Self(vec![])
    }

    /// The constant polynomial `1`.
    ///
    /// # Definition
    /// The multiplicative identity of `R[x]`, the image of `1` under the inclusion `R -> R[x]`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn one() -> Self {
        let one = R::default().one();
        Self(vec![one])
    }

    /// The degree of `self`, or `None` for zero polynomial.
    ///
    /// # Definition
    /// `deg f` is the greatest `k` with `f(k)` nonzero. The zero polynomial has no such `k`; the
    /// convention `deg 0 = -∞` is what makes `deg(fg) = deg f + deg g` hold over an integral
    /// domain, and `None` is that convention. By the invariant it is the stored length minus one.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn degree(&self) -> Option<usize> {
        let len = self.0.len();
        if len == 0 { None } else { Some(len - 1) }
    }

    /// Reduces `self` modulo `x^k`.
    ///
    /// # Definition
    /// Drops the coefficients of `x^k` and beyond, the `R`-linear projection of `R[x]` onto the
    /// span of `1, ..., x^{k-1}`. It is the quotient map `R[x] -> R[x]/(x^k)` followed by the
    /// representative of degree less than `k`, so it is additive but not multiplicative: the
    /// product of two reductions agrees with the reduction of the product only modulo `x^k`.
    ///
    /// # Complexity
    /// - Time: O(n - k)
    /// - Space: O(1)
    pub fn truncate(&mut self, k: usize) {
        self.0.truncate(k);
        self.normalize();
    }

    /// The formal derivative of `self`.
    ///
    /// # Definition
    /// `d/dx` is the `R`-linear map with `x^k -> k x^{k-1}`; it is a derivation,
    /// `(fg)' = f'g + fg'`, and uses no division, so it is defined over any semiring.
    /// The degree drops by one in characteristic zero, but may drop further otherwise:
    /// over `F_p` the derivative of `x^p` is zero.
    ///
    /// # Complexity
    /// - Time: O(n) multiplications in `R`
    /// - Space: O(n)
    pub fn derivative(&self) -> Self {
        let ring = R::default();
        let one = ring.one();
        let mut k = ring.zero();
        Self::from_vec(
            self.0
                .iter()
                .skip(1)
                .map(|v| {
                    k = ring.add(&k, &one);
                    ring.mul(&k, v)
                })
                .collect(),
        )
    }

    /// The stored coefficients `f(0), ..., f(deg)`, in increasing order.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn iter(&self) -> impl Iterator<Item = &R::Value> {
        self.0.iter()
    }

    /// The stored coefficients `f(0), ..., f(deg)`, in increasing order.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn into_vec(self) -> Vec<R::Value> {
        self.0
    }

    /// Restores the invariant by dropping the trailing zeros.
    ///
    /// # Complexity
    /// - Time: O(d), where `d` is the number of coefficients dropped
    /// - Space: O(1)
    fn normalize(&mut self) {
        let zero = R::default().zero();
        let len = self.0.iter().rposition(|v| *v != zero).map_or(0, |i| i + 1);
        self.0.truncate(len);
    }
}

impl<R: Semiring<Value: PartialEq + Clone> + Default> Poly<R> {
    /// Multiplies `self` by `x^j`.
    ///
    /// # Definition
    /// `f -> x^j f` raises every exponent by `j`. It mirrors `n << j = n 2^j` on the integers with
    /// `x` in place of the base: injective, `R`-linear, and not surjective, with `shr` as its left
    /// inverse.
    ///
    /// # Complexity
    /// - Time: O(n + j)
    /// - Space: O(n + j)
    pub fn shl(&self, j: usize) -> Self {
        if self.0.is_empty() {
            return Self::zero();
        }
        let ring = R::default();
        let mut value = Vec::with_capacity(self.0.len() + j);
        value.resize_with(j, || ring.zero());
        value.extend(self.0.iter().cloned());
        Self(value)
    }

    /// The quotient of `self` by `x^j`.
    ///
    /// # Definition
    /// `f -> f / x^j` drops the lowest `j` coefficients, mirroring `n >> j` on the integers. It is
    /// `R`-linear but neither injective nor multiplicative. `shl(j)` followed by `shr(j)` is the
    /// identity; the other order is the identity exactly when `x^j` divides `f`.
    ///
    /// # Complexity
    /// - Time: O(n - j)
    /// - Space: O(n - j)
    pub fn shr(&self, j: usize) -> Self {
        Self(self.0.get(j..).map_or_else(Vec::new, <[R::Value]>::to_vec))
    }

    /// The coefficients `f(0), ..., f(k - 1)`, padding with zero beyond the degree.
    pub fn coefficients(&self, k: usize) -> impl Iterator<Item = R::Value> + '_ {
        let ring = R::default();
        self.0
            .iter()
            .cloned()
            .chain(std::iter::repeat_with(move || ring.zero()))
            .take(k)
    }
}
impl<R: Field<Value: PartialEq> + Default> Poly<R> {
    /// The integral of `self` with zero constant term.
    ///
    /// # Definition
    /// The `R`-linear map with `x^k -> x^{k+1} / (k + 1)`, the section of `d/dx` picking the
    /// antiderivative that vanishes at `0`. It is a right inverse only: `(∫f)' = f` always, while
    /// `∫(f') = f - f(0)`, since the kernel of `d/dx` is the constants.
    ///
    /// # Contract
    /// `1, ..., deg + 1` are invertible in `R`, that is the characteristic is `0` or greater than
    /// `deg + 1`.
    ///
    /// # Complexity
    /// - Time: O(n) multiplications in `R` and one inversion
    /// - Space: O(n)
    pub fn integral(&self) -> Self {
        let ring = R::default();
        let n = self.0.len();
        if n == 0 {
            return Self::zero();
        }
        let one = ring.one();
        let mut ks = Vec::with_capacity(n);
        let mut fact = Vec::with_capacity(n);
        ks.push(ring.one());
        fact.push(ring.one());
        for i in 1..n {
            let k = ring.add(&ks[i - 1], &one);
            let f = ring.mul(&fact[i - 1], &k);
            ks.push(k);
            fact.push(f);
        }
        let mut value: Vec<R::Value> = (0..=n).map(|_| ring.zero()).collect();
        let mut acc = ring.inv(&fact[n - 1]);
        for i in (0..n).rev() {
            let inv_k = if i == 0 {
                ring.mul(&acc, &one)
            } else {
                ring.mul(&fact[i - 1], &acc)
            };
            value[i + 1] = ring.mul(&self.0[i], &inv_k);
            acc = ring.mul(&acc, &ks[i]);
        }
        Self(value)
    }
}

impl<R: RootOfUnity<Value: PartialEq + Clone> + Default> Poly<R> {
    /// The polynomial `f(x + c)`.
    ///
    /// # Definition
    /// The image of `f` under the `R`-algebra endomorphism of `R[x]` sending `x` to `x + c`. It is
    /// an automorphism, with inverse `x -> x - c`, so the degree and the leading coefficient are
    /// unchanged. Coefficientwise `[x^k] f(x+c) = Σ_{i >= k} binom(i, k) f(i) c^{i-k}`.
    ///
    /// It is not an operation on `R[[x]]`: the sum over `i >= k` is infinite there, and it does not
    /// descend to `R[x]/(x^k)` either, since `(x+c)^k` does not lie in `(x^k)`.
    ///
    /// # Contract
    /// `1, ..., deg` are invertible in `R`, that is the characteristic is `0` or greater than
    /// `deg`.
    ///
    /// # Complexity
    /// - Time: O(n log n)
    /// - Space: O(n)
    pub fn taylor_shift(&self, c: &R::Value) -> Self {
        let ring = R::default();
        let n = self.0.len();
        if n == 0 {
            return Self::zero();
        }
        let one = ring.one();
        let minus_one = ring.neg(&one);

        let mut g: Vec<R::Value> = Vec::with_capacity(n);
        let mut fact = ring.one();
        let mut last = ring.zero();
        for i in 0..n {
            if i > 0 {
                last = ring.add(&last, &one);
                fact = ring.mul(&fact, &last);
            }
            g.push(ring.mul(&fact, &self.0[i]));
        }
        g.reverse();
        let inv_fact = ring.inv(&fact);

        let mut h: Vec<R::Value> = Vec::with_capacity(n);
        let mut pow = ring.one();
        for _ in 0..n {
            h.push(ring.mul(&pow, &one));
            pow = ring.mul(&pow, c);
        }
        let mut acc = ring.mul(&inv_fact, &one);
        let mut k = ring.mul(&last, &one);
        for h in h.iter_mut().rev() {
            *h = ring.mul(h, &acc);
            acc = ring.mul(&acc, &k);
            k = ring.add(&k, &minus_one);
        }

        let mut p = Additive::default().convolve(&ring, g, h);
        p.truncate(n);
        let mut acc = ring.mul(&inv_fact, &one);
        let mut k = ring.mul(&last, &one);
        for p in p.iter_mut() {
            *p = ring.mul(p, &acc);
            acc = ring.mul(&acc, &k);
            k = ring.add(&k, &minus_one);
        }
        p.reverse();
        Self(p)
    }
}

impl<R: Semiring<Value: PartialEq + Clone> + Default> Clone for Poly<R> {
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    fn clone(&self) -> Self {
        Self(self.0.clone())
    }
}

impl<R: Ring<Value: PartialEq> + Default> std::ops::Neg for Poly<R> {
    type Output = Self;
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(1)
    fn neg(mut self) -> Self::Output {
        let ring = R::default();
        for v in self.0.iter_mut() {
            *v = ring.neg(v);
        }
        self
    }
}
impl<R: Ring<Value: PartialEq> + Default> std::ops::Neg for &Poly<R> {
    type Output = Poly<R>;
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    fn neg(self) -> Poly<R> {
        let ring = R::default();
        Poly(self.0.iter().map(|x| ring.neg(x)).collect())
    }
}

impl<R: Semiring<Value: PartialEq> + Default> std::ops::Add<&Self> for Poly<R> {
    type Output = Self;
    /// # Complexity
    /// - Time: O(n + m)
    /// - Space: O(1)
    fn add(mut self, rhs: &Self) -> Self {
        let ring = R::default();
        if self.0.len() < rhs.0.len() {
            self.0.resize_with(rhs.0.len(), || ring.zero());
        }
        for (l, r) in self.0.iter_mut().zip(rhs.0.iter()) {
            *l = ring.add(l, r);
        }
        self.normalize();
        self
    }
}
impl<R: Ring<Value: PartialEq> + Default> std::ops::Sub<&Self> for Poly<R> {
    type Output = Self;
    /// # Complexity
    /// - Time: O(n + m)
    /// - Space: O(1)
    fn sub(mut self, rhs: &Self) -> Self {
        let ring = R::default();
        if self.0.len() < rhs.0.len() {
            self.0.resize_with(rhs.0.len(), || ring.zero());
        }
        for (l, r) in self.0.iter_mut().zip(rhs.0.iter()) {
            *l = ring.add(l, &ring.neg(r));
        }
        self.normalize();
        self
    }
}
macro_rules! forward_binop {
    ($($trait:ident, $method:ident, $scalars:ident);* $(;)?) => {
        $(
            impl<R: $scalars<Value: PartialEq> + Default> std::ops::$trait for Poly<R> {
                type Output = Self;
                /// # Complexity
                /// - Time: O(n)
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

impl<R: Semiring<Value: PartialEq + Clone> + Default> std::ops::Add<Self> for &Poly<R> {
    type Output = Poly<R>;
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    fn add(self, rhs: Self) -> Poly<R> {
        let ring = R::default();
        let mut value = self.0.clone();
        if value.len() < rhs.0.len() {
            value.resize_with(rhs.0.len(), || ring.zero());
        }
        for (l, r) in value.iter_mut().zip(&rhs.0) {
            *l = ring.add(l, r);
        }
        let mut poly = Poly(value);
        poly.normalize();
        poly
    }
}
impl<R: Ring<Value: PartialEq + Clone> + Default> std::ops::Sub<Self> for &Poly<R> {
    type Output = Poly<R>;
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    fn sub(self, rhs: Self) -> Poly<R> {
        let ring = R::default();
        let mut value = self.0.clone();
        if value.len() < rhs.0.len() {
            value.resize_with(rhs.0.len(), || ring.zero());
        }
        for (l, r) in value.iter_mut().zip(&rhs.0) {
            *l = ring.add(l, &ring.neg(r));
        }
        let mut poly = Poly(value);
        poly.normalize();
        poly
    }
}
macro_rules! forward_ref_binop {
    ($($trait:ident, $method:ident, $scalars:ident);* $(;)?) => {
        $(
            impl<R: $scalars<Value: PartialEq + Clone> + Default> std::ops::$trait<Poly<R>> for &Poly<R> {
                type Output = Poly<R>;
                /// # Complexity
                /// - Time: O(n)
                /// - Space: O(n)
                fn $method(self, rhs: Poly<R>) -> Poly<R> {
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

impl<R: Ring<Value: PartialEq> + RootOfUnity + Default> std::ops::Mul<Self> for Poly<R> {
    type Output = Self;
    /// # Complexity
    /// - Time: O((n + m) log (n + m))
    /// - Space: O(n + m)
    fn mul(self, rhs: Self) -> Self {
        let ring = R::default();
        let h = Additive(Canonical::new()).convolve(&ring, self.0, rhs.0);
        Self(h)
    }
}
impl<R: Ring<Value: PartialEq + Clone> + RootOfUnity + Default> std::ops::Mul<&Self> for Poly<R> {
    type Output = Self;
    /// # Complexity
    /// - Time: O((n + m) log (n + m))
    /// - Space: O(n + m)
    fn mul(self, rhs: &Self) -> Self {
        self * rhs.clone()
    }
}
impl<R: Ring<Value: PartialEq + Clone> + RootOfUnity + Default> std::ops::Mul<Poly<R>>
    for &Poly<R>
{
    type Output = Poly<R>;
    /// # Complexity
    /// - Time: O((n + m) log (n + m))
    /// - Space: O(n + m)
    fn mul(self, rhs: Poly<R>) -> Poly<R> {
        self.clone() * rhs
    }
}
impl<R: Ring<Value: PartialEq + Clone> + RootOfUnity + Default> std::ops::Mul<Self> for &Poly<R> {
    type Output = Poly<R>;
    /// # Complexity
    /// - Time: O((n + m) log (n + m))
    /// - Space: O(n + m)
    fn mul(self, rhs: Self) -> Poly<R> {
        self.clone() * rhs.clone()
    }
}

impl<R: Semiring<Value: PartialEq> + Default> std::ops::AddAssign<&Self> for Poly<R> {
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(1)
    fn add_assign(&mut self, rhs: &Self) {
        let ring = R::default();
        if self.0.len() < rhs.0.len() {
            self.0.resize_with(rhs.0.len(), || ring.zero());
        }
        for (l, r) in self.0.iter_mut().zip(rhs.0.iter()) {
            *l = ring.add(l, r);
        }
        self.normalize();
    }
}
impl<R: Ring<Value: PartialEq> + Default> std::ops::SubAssign<&Self> for Poly<R> {
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(1)
    fn sub_assign(&mut self, rhs: &Self) {
        let ring = R::default();
        if self.0.len() < rhs.0.len() {
            self.0.resize_with(rhs.0.len(), || ring.zero());
        }
        for (l, r) in self.0.iter_mut().zip(rhs.0.iter()) {
            *l = ring.add(l, &ring.neg(r));
        }
        self.normalize();
    }
}
macro_rules! forward_op_assign {
    ($($trait:ident, $method:ident, $scalars:ident);* $(;)?) => {
        $(
            impl<R: $scalars<Value: PartialEq> + Default> std::ops::$trait for Poly<R> {
                /// # Complexity
                /// - Time: O(n)
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

impl<R: Ring<Value: PartialEq> + RootOfUnity + Default> std::ops::MulAssign<Self> for Poly<R> {
    /// # Complexity
    /// - Time: O((n + m) log (n + m))
    /// - Space: O(n + m)
    fn mul_assign(&mut self, rhs: Self) {
        let lhs = std::mem::take(&mut self.0);
        self.0 = Additive::default().convolve(&R::default(), lhs, rhs.0);
    }
}
impl<R: Ring<Value: PartialEq + Clone> + RootOfUnity + Default> std::ops::MulAssign<&Self>
    for Poly<R>
{
    /// # Complexity
    /// - Time: O((n + m) log (n + m))
    /// - Space: O(n + m)
    fn mul_assign(&mut self, rhs: &Self) {
        *self *= rhs.clone();
    }
}

impl<R: Semiring<Value: PartialEq> + Default> std::ops::Index<usize> for Poly<R> {
    type Output = R::Value;
    /// # Panics
    /// Panics if `index` is out of bounds.
    fn index(&self, index: usize) -> &Self::Output {
        assert!(
            index < self.0.len(),
            "index out of bounds: index={index}, len={}",
            self.0.len()
        );
        &self.0[index]
    }
}
