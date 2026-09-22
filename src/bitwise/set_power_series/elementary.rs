use crate::algebra::Ring;
use crate::bitwise::set_power_series::{ranked_subset_mobius, ranked_subset_zeta, sps_convolve};

/// The exponential of `f`, as a set power series.
///
/// # Definition
/// `(exp f)(S) = Σ_π Π_{B∈π} f(B)`, the sum over the partitions `π` of `S` into nonempty blocks,
/// with `(exp f)(∅) = 1`. If `1, ..., n` are invertible in `R`, it is `Σ_{k=0,...,n} f^k / k!`.
///
/// # Contract
/// `f(∅) = 0`.
///
/// # Complexity
/// - Time: O(2^n n^2)
/// - Space: O(2^n n)
///
/// # Panics
/// Panics if `f.len()` is not a power of two.
pub fn sps_exp<R: Ring<Value: Clone>>(ring: &R, f: &[R::Value]) -> Vec<R::Value> {
    assert!(
        f.len().is_power_of_two(),
        "length must be power of two: len={}",
        f.len()
    );
    let n = f.len().trailing_zeros() as usize;
    let mut g = Vec::with_capacity(f.len());
    g.push(ring.one());
    for i in 0..n {
        let w = 1 << i;
        let high = sps_convolve(ring, f[w..w << 1].to_vec(), g[..w].to_vec());
        g.extend(high);
    }
    g
}

/// The logarithm of `f`, as a set power series, inverting [`sps_exp`].
///
/// # Definition
/// The unique `g` with `g(∅) = 0` and `exp g = f`, which exists as `exp` is a bijection from
/// `{g: g(∅) = 0}` onto `{f: f(∅) = 1}`.
///
/// # Contract
/// `f(∅) = 1`.
///
/// # Complexity
/// - Time: O(2^n n^2)
/// - Space: O(2^n n)
///
/// # Panics
/// Panics if `f.len()` is not a power of two.
pub fn sps_log<R: Ring<Value: Clone>>(ring: &R, f: &[R::Value]) -> Vec<R::Value> {
    assert!(
        f.len().is_power_of_two(),
        "length must be power of two: len={}",
        f.len(),
    );
    let n = f.len().trailing_zeros() as usize;
    let mut g = Vec::with_capacity(f.len());
    g.push(ring.zero());
    for i in 0..n {
        let w = 1 << i;
        let mut a = ranked_subset_zeta(ring, f[w..w << 1].to_vec());
        let b = ranked_subset_zeta(ring, f[..w].to_vec());
        for s in 0..w {
            for d in 0..=i {
                for k in 0..d {
                    let x = ring.mul(&a[k][s], &b[d - k][s]);
                    a[d][s] = ring.add(&a[d][s], &ring.neg(&x));
                }
            }
        }
        g.extend(ranked_subset_mobius(ring, a));
    }
    g
}
