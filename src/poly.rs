use crate::algebra::{Field, Ring, RootOfUnity, Semiring};
use crate::periodic_function::cyclic_convolve;

/// The polynomial ring `R[x]` over a semiring `R`.
///
/// # Definition
/// The free `R`-module on the monomials `x^k`, `k >= 0`, with the product
/// `x^i x^j = x^{i+j}` extended bilinearly: `(fg)(k) = Σ_{i+j=k} f(i) g(j)`. It is the monoid
/// algebra `R[(N, +)]`.
///
/// # Invariants
/// - The last stored coefficient is nonzero, so the zero polynomial is the empty vector and the
///   stored length is `deg + 1`.
pub struct Poly<R: Semiring<Value: PartialEq> + Default>(Vec<R::Value>);

impl<R: Semiring<Value: PartialEq> + Default> Poly<R> {
    /// Constructs `Σ_k value[k] x^k`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(1)
    pub fn from_vec(value: Vec<R::Value>) -> Self {
        let mut poly = Self(value);
        poly.normalize();
        poly
    }

    /// The zero polynomial, the additive identity of `R[x]`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn zero() -> Self {
        Self(vec![])
    }

    /// The constant polynomial `1`, the multiplicative identity of `R[x]`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn one() -> Self {
        Self(vec![R::default().one()])
    }

    /// The degree of `self`.
    ///
    /// # Definition
    /// The greatest `k` with `f(k) != 0`, or `None` for the zero polynomial.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn degree(&self) -> Option<usize> {
        self.0.len().checked_sub(1)
    }

    /// Reduces `self` modulo `x^k`, in place.
    ///
    /// # Definition
    /// `f -> Σ_{i<k} f(i) x^i`, the representative of degree less than `k` of `f` in
    /// `R[x]/(x^k)`.
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
    /// The `R`-linear map `x^k -> k x^{k-1}`. Over `F_p` it may lower the degree by more than one,
    /// since `k = 0` for `p | k`.
    ///
    /// # Complexity
    /// - Time: O(n)
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

    /// The coefficients `f(0), ..., f(deg)`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn iter(&self) -> impl Iterator<Item = &R::Value> {
        self.0.iter()
    }

    /// The coefficients `f(0), ..., f(deg)`.
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
    /// The product `x^j f`.
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
    /// `f -> Σ_{i>=j} f(i) x^{i-j}`, dropping the coefficients of degree less than `j`.
    ///
    /// # Complexity
    /// - Time: O(n - j)
    /// - Space: O(n - j)
    pub fn shr(&self, j: usize) -> Self {
        Self(self.0.get(j..).map_or_else(Vec::new, <[R::Value]>::to_vec))
    }

    /// The coefficients `f(0), ..., f(k - 1)`, with `f(i) = 0` for `i > deg`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
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
    /// The integral of `self` vanishing at `0`.
    ///
    /// # Definition
    /// The `R`-linear map `x^k -> x^{k+1} / (k + 1)`, the right inverse of the derivative with zero
    /// constant term.
    ///
    /// # Contract
    /// `1, ..., deg + 1` are invertible in `R`, that is the characteristic is `0` or greater than
    /// `deg + 1`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    pub fn integral(&self) -> Self {
        let ring = R::default();
        let n = self.0.len();
        if n == 0 {
            return Self::zero();
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
            value[i + 1] = ring.mul(&self.0[i], &ring.mul(&fact[i], &acc));
            acc = ring.mul(&acc, &k);
            k = ring.add(&k, &minus_one);
        }
        Self(value)
    }
}

impl<R: RootOfUnity<Value: PartialEq> + Default> Poly<R> {
    /// The polynomial `f(x + c)`.
    ///
    /// # Definition
    /// The image of `f` under the `R`-algebra automorphism of `R[x]` sending `x` to `x + c`:
    /// `[x^k] f(x + c) = Σ_{i>=k} binom(i, k) f(i) c^{i-k}`.
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
        let mut k = ring.zero();
        for i in 0..n {
            if i > 0 {
                k = ring.add(&k, &one);
                fact = ring.mul(&fact, &k);
            }
            g.push(ring.mul(&fact, &self.0[i]));
        }
        g.reverse();

        let mut inv_fact: Vec<R::Value> = Vec::with_capacity(n);
        inv_fact.push(ring.inv(&fact));
        for _ in 1..n {
            inv_fact.push(ring.mul(inv_fact.last().unwrap(), &k));
            k = ring.add(&k, &minus_one);
        }
        inv_fact.reverse();

        let mut h: Vec<R::Value> = Vec::with_capacity(n);
        h.push(one);
        for j in 1..n {
            h.push(ring.mul(&h[j - 1], c));
        }
        for (h, inv) in h.iter_mut().zip(&inv_fact) {
            *h = ring.mul(h, inv);
        }

        let mut p = poly_convolve(&ring, g, h);
        p.truncate(n);
        for (p, inv) in p.iter_mut().zip(inv_fact.iter().rev()) {
            *p = ring.mul(p, inv);
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

impl<R: Semiring<Value: PartialEq> + Default> std::ops::AddAssign<&Self> for Poly<R> {
    /// # Complexity
    /// - Time: O(n + m)
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
    /// - Time: O(n + m)
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

impl<R: Ring<Value: PartialEq> + RootOfUnity + Default> std::ops::MulAssign<Self> for Poly<R> {
    /// # Complexity
    /// - Time: O((n + m) log (n + m))
    /// - Space: O(n + m)
    fn mul_assign(&mut self, rhs: Self) {
        let lhs = std::mem::take(&mut self.0);
        self.0 = poly_convolve(&R::default(), lhs, rhs.0);
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

macro_rules! forward_binop {
    ($($trait:ident, $method:ident, $assign:ident, $assign_method:ident, $scalars:ident);* $(;)?) => {
        $(
            impl<R: $scalars<Value: PartialEq> + Default> std::ops::$trait for Poly<R> {
                type Output = Self;
                fn $method(mut self, rhs: Self) -> Self {
                    std::ops::$assign::$assign_method(&mut self, rhs);
                    self
                }
            }
            impl<R: $scalars<Value: PartialEq + Clone> + Default> std::ops::$trait<&Self> for Poly<R> {
                type Output = Self;
                fn $method(mut self, rhs: &Self) -> Self {
                    std::ops::$assign::$assign_method(&mut self, rhs);
                    self
                }
            }
            impl<R: $scalars<Value: PartialEq + Clone> + Default> std::ops::$trait<Poly<R>> for &Poly<R> {
                type Output = Poly<R>;
                fn $method(self, rhs: Poly<R>) -> Poly<R> {
                    std::ops::$trait::$method(self.clone(), rhs)
                }
            }
            impl<R: $scalars<Value: PartialEq + Clone> + Default> std::ops::$trait<Self> for &Poly<R> {
                type Output = Poly<R>;
                fn $method(self, rhs: Self) -> Poly<R> {
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

/// The length of the shorter operand up to which [`poly_convolve`] multiplies directly.
const NAIVE_LIMIT: usize = 32;
/// The product of `f` and `g` in `R[x]`, of length `f.len() + g.len() - 1`.
///
/// # Definition
/// `(fg)(k) = Σ_{i+j=k} f(i) g(j)`, the product of the monoid algebra `R[(N, +)] = R[x]`. The
/// quotient map `R[x] -> R[x]/(x^n - 1) = R[Z/nZ]` is injective in degree less than `n`, so for
/// `n > deg fg` the product is that of `R[Z/nZ]` (see [`cyclic_convolve`]).
///
/// # Complexity
/// - Time: O((n + m) log (n + m))
/// - Space: O(n + m)
///
/// # Panics
/// Panics if `R` has no primitive `n`-th root of unity for the least power of two at least
/// `f.len() + g.len() - 1`.
pub fn poly_convolve<R: RootOfUnity>(
    ring: &R,
    mut f: Vec<R::Value>,
    mut g: Vec<R::Value>,
) -> Vec<R::Value> {
    if f.is_empty() || g.is_empty() {
        return Vec::new();
    }
    let len = f.len() + g.len() - 1;
    if f.len().min(g.len()) <= NAIVE_LIMIT {
        let mut h: Vec<R::Value> = (0..len).map(|_| ring.zero()).collect();
        for (i, fi) in f.iter().enumerate() {
            for (hj, gj) in h[i..].iter_mut().zip(g.iter()) {
                *hj = ring.add(hj, &ring.mul(fi, gj));
            }
        }
        return h;
    }
    let n = len.next_power_of_two();
    f.resize_with(n, || ring.zero());
    g.resize_with(n, || ring.zero());
    let mut h = cyclic_convolve(ring, f, g);
    h.truncate(len);
    h
}
