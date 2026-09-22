use crate::algebra::canonical::Canonical;
use crate::algebra::{Commutative, Inv, One, Zero};
use crate::arithmetic::prime::is_mersenne_exponent;

/// An element of the prime field `Z/MZ` for a Mersenne prime `M = 2^K - 1`, whose products are
/// reduced without division.
///
/// # Definition
/// `(Z/MZ, 0, 1, +, *)` is the prime field of order `M`. Identifying `[0, M)` with the residues,
/// the operations are `(a + b) % M` and `(a * b) % M`. `M` is prime and `K <= 63`, that is `K` is
/// one of `2, 3, 5, 7, 13, 17, 19, 31, 61`; any other `K` is rejected at compile time.
///
/// # Invariants
/// - `self.0` is in `[0, M)`.
///
/// # Complexity
/// - Space: O(1)
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Mersenne<const K: u32>(u64);

impl<const K: u32> Mersenne<K> {
    /// The modulus `M`, rejected at compile time unless it is prime and `K <= 63`.
    const MODULUS: u64 = {
        assert!(K <= 63, "K must be at most 63");
        assert!(is_mersenne_exponent(K), "2^K - 1 must be prime");
        (1 << K) - 1
    };

    /// The residue class of `n`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn new(n: u64) -> Self {
        Self(n % Self::MODULUS)
    }

    /// The residue of `x` modulo `M`, in `[0, M)`.
    ///
    /// # Contract
    /// `x <= (M - 1)^2`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    fn reduce(x: u128) -> u64 {
        let y = (x >> K) as u64 + (x as u64 & Self::MODULUS);
        if y >= Self::MODULUS {
            y - Self::MODULUS
        } else {
            y
        }
    }

    /// The residue of `self` in `[0, M)`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn val(self) -> u64 {
        self.0
    }

    /// The power `self^exp`.
    ///
    /// # Definition
    /// `self^0 = 1` and `self^e = self^{e-1} self` for `e >= 1`.
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

    /// The inverse `self^{-1}`.
    ///
    /// # Complexity
    /// - Time: O(K)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `self == 0`.
    pub fn inv(self) -> Self {
        assert!(self.0 != 0, "zero has no inverse in Z/(2^{K}-1)Z");
        self.pow(Self::MODULUS - 2)
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

impl<const K: u32> Commutative for Canonical<Mersenne<K>> {}

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
