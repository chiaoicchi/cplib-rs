use crate::algebra::canonical::Canonical;
use crate::algebra::{Commutative, Inv, One, RootOfUnity, Zero};
use crate::arithmetic::prime::{is_prime, primitive_root};

/// An element of the prime field `Z/PZ`.
///
/// # Definition
/// `(Z/PZ, 0, 1, +, *)` is the prime field of order `P`. Identifying `[0, P)` with the residues,
/// the operations are `(a + b) % P` and `(a * b) % P`. `P` is a prime less than `2^31`; any other
/// `P` is rejected at compile time.
///
/// # Invariants
/// - `self.0` is in `[0, P)`.
///
/// # Complexity
/// - Space: O(1)
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Fp<const P: u32>(u32);

impl<const P: u32> Fp<P> {
    /// The modulus `P`, rejected at compile time unless it is a prime less than `2^31`.
    const MODULUS: u32 = {
        assert!(is_prime(P), "P must be prime");
        assert!(P < 1 << 31, "P must be less than 2^31");
        P
    };

    /// The least generator of the multiplicative group `(Z/PZ)^*`, as its residue.
    pub const PRIMITIVE_ROOT: u32 = primitive_root(Self::MODULUS);

    /// The residue class of `n`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn new(n: u32) -> Self {
        Self(n % Self::MODULUS)
    }

    /// The residue of `self` in `[0, P)`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn val(self) -> u32 {
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
        let mut x = Fp(1);
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
    /// - Time: O(log P)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `self == 0`.
    pub fn inv(self) -> Self {
        assert!(self.0 != 0, "zero has no inverse in Z/{P}Z");
        self.pow(P as u64 - 2)
    }
}

#[doc(hidden)]
#[macro_export]
macro_rules! __fp {
    ($value:expr) => {
        $crate::num::fp::Fp::from($value)
    };
    ($value:expr, mod $p:expr) => {
        $crate::num::fp::Fp::<{ $p }>::from($value)
    };
}
/// The element of [`Fp`] given by `x`.
///
/// # Definition
/// `fp!(x)` is `Fp::from(x)`, and `fp!(x, mod p)` is `Fp::<p>::from(x)`.
pub use __fp as fp;
macro_rules! impl_from_signed {
    ($($t:ty => $w:ty), * $(,)?) => {$(
        impl<const P: u32> From<$t> for Fp<P> {
            fn from(x: $t) -> Self {
                Self((x as $w).rem_euclid(Self::MODULUS as $w) as u32)
            }
        }
    )*};
}
impl_from_signed!(i8 => i32, i16 => i32, i32 => i32, i64 => i64, i128 => i128, isize => i64);
macro_rules! impl_from_unsigned {
    ($($t:ty => $w:ty), *) => {$(
        impl<const P: u32> From<$t> for Fp<P> {
            fn from(x: $t) -> Self {
                Self((x as $w % Self::MODULUS as $w) as u32)
            }
        }
    )*}
}
impl_from_unsigned!(u8 => u32, u16 => u32, u32 => u32, u64 => u64, u128 => u128, usize => u64);

impl<const P: u32> std::ops::Neg for Fp<P> {
    type Output = Self;
    #[inline]
    fn neg(mut self) -> Self::Output {
        if self.0 > 0 {
            self.0 = P - self.0;
        }
        self
    }
}
impl<const P: u32> std::ops::Neg for &Fp<P> {
    type Output = Fp<P>;
    #[inline]
    fn neg(self) -> Fp<P> {
        -*self
    }
}

impl<const P: u32> std::ops::Add for Fp<P> {
    type Output = Self;
    #[inline]
    fn add(mut self, rhs: Self) -> Self {
        self.0 += rhs.0;
        if self.0 >= P {
            self.0 -= P;
        }
        self
    }
}
impl<const P: u32> std::ops::Sub for Fp<P> {
    type Output = Self;
    #[inline]
    fn sub(mut self, rhs: Self) -> Self {
        if self.0 < rhs.0 {
            self.0 += P;
        }
        self.0 -= rhs.0;
        self
    }
}
impl<const P: u32> std::ops::Mul for Fp<P> {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self((self.0 as u64 * rhs.0 as u64 % P as u64) as u32)
    }
}
#[allow(clippy::suspicious_arithmetic_impl)]
impl<const P: u32> std::ops::Div for Fp<P> {
    type Output = Self;
    #[inline]
    fn div(self, rhs: Self) -> Self {
        self * rhs.inv()
    }
}
macro_rules! forward_ref_binop {
    ($($trait:ident, $method:ident);* $(;)?) => {
        $(
            impl<const P:u32> std::ops::$trait<&Self> for Fp<P> {
                type Output = Self;
                #[inline]
                fn $method(self, rhs: &Self) -> Self {
                    std::ops::$trait::$method(self, *rhs)
                }
            }

            impl<const P:u32> std::ops::$trait<Fp<P>> for &Fp<P> {
                type Output = Fp<P>;
                #[inline]
                fn $method(self, rhs: Fp<P>) -> Fp<P> {
                    std::ops::$trait::$method(*self, rhs)
                }
            }

            impl<const P: u32> std::ops::$trait for &Fp<P> {
                type Output = Fp<P>;
                #[inline]
                fn $method(self, rhs: Self) -> Fp<P> {
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

impl<const P: u32> std::ops::AddAssign for Fp<P> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0;
        if self.0 >= P {
            self.0 -= P;
        }
    }
}
impl<const P: u32> std::ops::SubAssign for Fp<P> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        if self.0 < rhs.0 {
            self.0 += P;
        }
        self.0 -= rhs.0;
    }
}
impl<const P: u32> std::ops::MulAssign for Fp<P> {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        self.0 = (self.0 as u64 * rhs.0 as u64 % P as u64) as u32;
    }
}
impl<const P: u32> std::ops::DivAssign for Fp<P> {
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}
macro_rules! forward_ref_op_assign {
    ($($trait:ident, $method:ident);* $(;)?) => {
        $(
            impl<const P:u32> std::ops::$trait<&Self> for Fp<P> {
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
        impl<const P: u32> std::iter::$trait for Fp<P> {
            fn $method<I: Iterator<Item = Self>>(iter: I) -> Self {
                iter.fold($id, |acc, x| acc $op x)
            }
        }
        impl<'a, const P: u32> std::iter::$trait<&'a Fp<P>> for Fp<P> {
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

impl<const P: u32> Zero for Fp<P> {
    fn zero() -> Self {
        Self::new(0)
    }
}
impl<const P: u32> One for Fp<P> {
    fn one() -> Self {
        Self::new(1)
    }
}
impl<const P: u32> Inv for Fp<P> {
    type Output = Fp<P>;
    fn inv(&self) -> Self {
        Self::inv(*self)
    }
}
impl<const P: u32> Commutative for Canonical<Fp<P>> {}
impl<const P: u32> RootOfUnity for Canonical<Fp<P>> {
    fn root_of_unity(&self, n: usize) -> Option<Fp<P>> {
        let m = (P - 1) as usize;
        if n == 0 || m % n != 0 {
            return None;
        }
        Some(Fp::new(Fp::<P>::PRIMITIVE_ROOT).pow((m / n) as u64))
    }
}

impl<const P: u32> std::fmt::Debug for Fp<P> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}
impl<const P: u32> std::fmt::Display for Fp<P> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}
