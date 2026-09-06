use crate::algebra::{Inv, Zero};

/// The greatest common divisor by the Euclidean algorithm.
///
/// # Definition
/// `gcd(a, 0) = a` and `gcd(a, b) = gcd(b, a mod b)` for `b != 0`; in particular
/// `gcd(0, 0) = 0`.
///
/// # Complexity
/// - Time: O(log min(a, b))
/// - Space: O(1)
pub fn gcd<T: Clone + PartialEq + Zero + std::ops::Rem<Output = T>>(a: T, b: T) -> T {
    let (mut a, mut b) = (a, b);
    while b != T::zero() {
        (a, b) = (b.clone(), a % b);
    }
    a
}

/// The least common multiple by the Euclidean algorithm.
///
/// # Definition
/// `lcm(0, 0) = 0` and `lcm(a, b) * gcd(a, b) = a * b`.
///
/// # Complexity
/// - Time: O(log min(a, b))
/// - Space: O(1)
pub fn lcm<
    T: Clone
        + PartialEq
        + Zero
        + std::ops::Mul<Output = T>
        + std::ops::Div<Output = T>
        + std::ops::Rem<Output = T>,
>(
    a: T,
    b: T,
) -> T {
    a.clone() / gcd(a, b.clone()) * b
}
