use crate::algebra::{Inv, One, Zero};
use crate::num::prime::is_prime;

/// An element of the prime field `Z/(2^K - 1)Z` for a Mersenne prime `2^K - 1`.
///
/// # Definition
/// `Z/(2^K - 1)Z` is the quotient ring of the integers modulo the Mersenne number `2^K - 1`,
/// with addition and multiplication induced from `Z`. When `2^K - 1` is prime it is a field,
/// and `2^K = 1 (mod 2^K - 1)` lets a product `x < 2^(2K)` be reduced without division, as
/// `x = a 2^K + b = a + b (mod 2^K - 1)`.
///
/// # Invariants
/// - `2^K - 1` is prime and `K <= 63`, so that the product of two residues fits in `u128`.
///   A `K` violating this is rejected at compile time.
/// - An element is represented by its canonical residue in `[0, 2^K - 1)`.
///
/// # Complexity
/// - Space: O(1)
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Mersenne<const K: u32>(u64);

impl<const K: u32> Mersenne<K> {
    /// The modulus `2^K - 1`, checked at compile time to be prime with `K <= 63`.
    const MODULUS: u64 = {
        assert!(K <= 63, "K must be at most 63");
        assert!(is_mersenne_prime(K), "2^K - 1 must be prime");
        (1 << K) - 1
    };

    /// Creates an element in Mersenne from a value `n` reduced modulo `2^K-1`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn new(n: u64) -> Self {
        Self(Self::reduce(n as u128))
    }

    /// Reduces `x < 2^{2K}` modulo `2^K - 1`, using `2^K = 1 (mod 2^K - 1)`:
    /// `x = a 2^K + b = a + b (mod 2^K - 1)`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    fn reduce(x: u128) -> u64 {
        let m = Self::MODULUS as u128;
        let y = (x >> K) + (x & m);
        let y = (y >> K) + (y & m);
        let y = y as u64;
        if y >= Self::MODULUS {
            y - Self::MODULUS
        } else {
            y
        }
    }

    /// Returns the inner value in `[0, 2^K - 1)`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn value(self) -> u64 {
        self.0
    }

    /// Raises `self` to the power of `exp`.
    ///
    /// # Complexity
    /// - Time: O(log exp)
    /// - Space: O(1)
    pub fn pow(mut self, mut exp: u64) -> Self {
        let mut x = Mersenne(1);
        while exp > 0 {
            if exp & 1 == 1 {
                x *= self;
            }
            self *= self;
            exp >>= 1;
        }
        x
    }

    /// Returns the multiplicative inverse of `self`,
    /// i.e. the element `x` with `self * x = x * self = 1`.
    ///
    /// # Complexity
    /// - Time: O(K)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `self == 0`.
    pub fn inv(self) -> Self {
        assert!(self.0 != 0, "zero has no inverse in Z/(2^{K}-1)Z");
        self.pow(Self::MODULUS as u64 - 2)
    }
}

impl<const K: u32> std::ops::Neg for Mersenne<K> {
    type Output = Self;
    #[inline]
    fn neg(mut self) -> Self::Output {
        if self.0 > 0 {
            self.0 = Self::MODULUS - self.0;
        }
        self
    }
}
impl<const K: u32> std::ops::Neg for &Mersenne<K> {
    type Output = Mersenne<K>;
    #[inline]
    fn neg(self) -> Mersenne<K> {
        -*self
    }
}

impl<const K: u32> std::ops::Add for Mersenne<K> {
    type Output = Self;
    #[inline]
    fn add(mut self, rhs: Self) -> Self {
        self.0 += rhs.0;
        if self.0 >= Self::MODULUS {
            self.0 -= Self::MODULUS;
        }
        self
    }
}
impl<const K: u32> std::ops::Sub for Mersenne<K> {
    type Output = Self;
    #[inline]
    fn sub(mut self, rhs: Self) -> Self {
        if self.0 < rhs.0 {
            self.0 += Self::MODULUS;
        }
        self.0 -= rhs.0;
        self
    }
}
impl<const K: u32> std::ops::Mul for Mersenne<K> {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self(Self::reduce(self.0 as u128 * rhs.0 as u128))
    }
}
#[allow(clippy::suspicious_arithmetic_impl)]
impl<const K: u32> std::ops::Div for Mersenne<K> {
    type Output = Self;
    #[inline]
    fn div(self, rhs: Self) -> Self {
        self * rhs.inv()
    }
}
macro_rules! forward_ref_binop {
    ($($trait:ident, $method:ident);* $(;)?) => {
        $(
            impl<const K:u32> std::ops::$trait<&Self> for Mersenne<K> {
                type Output = Self;
                #[inline]
                fn $method(self, rhs: &Self) -> Self {
                    std::ops::$trait::$method(self, *rhs)
                }
            }

            impl<const K:u32> std::ops::$trait<Mersenne<K>> for &Mersenne<K> {
                type Output = Mersenne<K>;
                #[inline]
                fn $method(self, rhs: Mersenne<K>) -> Mersenne<K> {
                    std::ops::$trait::$method(*self, rhs)
                }
            }

            impl<const K: u32> std::ops::$trait for &Mersenne<K> {
                type Output = Mersenne<K>;
                #[inline]
                fn $method(self, rhs: Self) -> Mersenne<K> {
                    std::ops::$trait::$method(*self, *rhs)
                }
            }
        )*
    };
}
forward_ref_binop! {
    Add, add;
    Sub, sub;
    Mul, mul;
    Div, div;
}

impl<const K: u32> std::ops::AddAssign for Mersenne<K> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0;
        if self.0 >= Self::MODULUS {
            self.0 -= Self::MODULUS;
        }
    }
}
impl<const K: u32> std::ops::SubAssign for Mersenne<K> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        if self.0 < rhs.0 {
            self.0 += Self::MODULUS;
        }
        self.0 -= rhs.0;
    }
}
impl<const K: u32> std::ops::MulAssign for Mersenne<K> {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        self.0 = Self::reduce(self.0 as u128 * rhs.0 as u128);
    }
}
impl<const K: u32> std::ops::DivAssign for Mersenne<K> {
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}
macro_rules! forward_ref_op_assign {
    ($($trait:ident, $method:ident);* $(;)?) => {
        $(
            impl<const K:u32> std::ops::$trait<&Self> for Mersenne<K> {
                #[inline]
                fn $method(&mut self, rhs: &Self) {
                    std::ops::$trait::$method(self, *rhs);
                }
            }
        )*
    };
}
forward_ref_op_assign! {
    AddAssign, add_assign;
    SubAssign, sub_assign;
    MulAssign, mul_assign;
    DivAssign, div_assign;
}

macro_rules! impl_fold {
    ($($trait:ident, $method:ident, $id:expr, $op:tt);* $(;)?) => {$(
        impl<const K: u32> std::iter::$trait for Mersenne<K> {
            fn $method<I: Iterator<Item = Self>>(iter: I) -> Self {
                iter.fold($id, |acc, x| acc $op x)
            }
        }
        impl<'a, const K: u32> std::iter::$trait<&'a Mersenne<K>> for Mersenne<K> {
            fn $method<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
                iter.fold($id, |acc, x| acc $op x)
            }
        }
    )*};
}
impl_fold! {
    Sum, sum, Self::zero(), +;
    Product, product, Self::one(), *;
}

impl<const K: u32> Zero for Mersenne<K> {
    fn zero() -> Self {
        Self::new(0)
    }
}
impl<const K: u32> One for Mersenne<K> {
    fn one() -> Self {
        Self::new(1)
    }
}
impl<const K: u32> Inv for Mersenne<K> {
    type Output = Mersenne<K>;
    fn inv(&self) -> Self {
        Self::inv(*self)
    }
}

impl<const K: u32> std::fmt::Debug for Mersenne<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}
impl<const K: u32> std::fmt::Display for Mersenne<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// Returns `true` if `2^K - 1` is prime.
const fn is_mersenne_prime(k: u32) -> bool {
    if k == 2 {
        return true;
    }
    if k < 2 || !is_prime(k) {
        return false;
    }
    let m = (1u128 << k) - 1;
    let mut s = 4;
    let mut i = 0;
    while i < k - 2 {
        s = (s * s + m - 2) % m;
        i += 1;
    }
    s == 0
}
