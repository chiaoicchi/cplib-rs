use crate::algebra::{Bounded, One, Zero};

macro_rules! impl_zero_one {
    ($($t:ty),* $(,)?) => {$(
        impl Zero for $t {
            fn zero() -> Self {
                0
            }
        }
        impl One for $t {
            fn one() -> Self {
                1
            }
        }
    )*};
}
impl_zero_one!(
    i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize
);

impl<T: Zero> Zero for std::num::Wrapping<T> {
    fn zero() -> Self {
        Self(T::zero())
    }
}
impl<T: One> One for std::num::Wrapping<T> {
    fn one() -> Self {
        Self(T::one())
    }
}

macro_rules! impl_bounded {
    ($($t:ty),* $(,)?) => {$(
        impl Bounded for $t {
            fn min_value() -> Self {
                <$t>::MIN
            }
            fn max_value() -> Self {
                <$t>::MAX
            }
        }
    )*};
}
impl_bounded!(
    i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize
);

impl<T: Bounded> Bounded for std::num::Wrapping<T> {
    fn min_value() -> Self {
        Self(T::min_value())
    }
    fn max_value() -> Self {
        Self(T::max_value())
    }
}
