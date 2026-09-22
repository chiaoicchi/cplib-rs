pub mod elementary;

use crate::algebra::Ring;
use crate::bitwise::or::{subset_mobius, subset_zeta};

/// The product of `f` and `g` as set power series, the subset convolution.
///
/// # Definition
/// `(fg)(S) = Σ_{T⊆S} f(T) g(S\T)`.
///
/// # Complexity
/// - Time: O(2^n n^2)
/// - Space: O(2^n n)
///
/// # Panics
/// Panics if `f` and `g` differ in length.
/// Panics if that length is not a power of two.
pub fn sps_convolve<R: Ring>(ring: &R, f: Vec<R::Value>, g: Vec<R::Value>) -> Vec<R::Value> {
    assert_eq!(
        f.len(),
        g.len(),
        "f and g must have the same length: f={}, g={}",
        f.len(),
        g.len()
    );
    assert!(
        f.len().is_power_of_two(),
        "length must be a power of two: len={}",
        f.len()
    );
    let n = f.len().trailing_zeros() as usize;
    let (mut f, g) = (ranked_subset_zeta(ring, f), ranked_subset_zeta(ring, g));
    for s in 0..1 << n {
        for d in (0..=n).rev() {
            let mut acc = ring.mul(&f[d][s], &g[0][s]);
            for i in (0..d).rev() {
                acc = ring.add(&acc, &ring.mul(&f[i][s], &g[d - i][s]));
            }
            f[d][s] = acc;
        }
    }
    ranked_subset_mobius(ring, f)
}

/// The subset zeta transform of `f`, split by rank.
///
/// # Definition
/// The `n + 1` set functions `F_i(S) = Σ_{T⊆S, |T|=i} f(T)`, for `i` in `[0, n]`.
///
/// # Complexity
/// - Time: O(2^n n^2)
/// - Space: O(2^n n)
///
/// # Panics
/// Panics if `f.len()` is not a power of two.
pub fn ranked_subset_zeta<R: Ring>(ring: &R, f: Vec<R::Value>) -> Vec<Vec<R::Value>> {
    assert!(
        f.len().is_power_of_two(),
        "length must be a power of two: len={}",
        f.len()
    );
    let n = f.len().trailing_zeros();
    let mut layers: Vec<Vec<R::Value>> = (0..=n)
        .map(|_| (0..f.len()).map(|_| ring.zero()).collect())
        .collect();
    for (s, x) in f.into_iter().enumerate() {
        layers[s.count_ones() as usize][s] = x;
    }
    for layer in &mut layers {
        subset_zeta(ring, layer);
    }
    layers
}

/// The set function recoverd from its ranked zeta transform, inverting [`ranked_subset_zeta`].
///
/// # Definition
/// `f(S) = Σ_{T⊆S} (-1)^{|S|-|T|} F_{|S|}(T)`, for the set functions `F_0, ..., F_n` given as
/// `layers`.
///
/// # Contract
/// `layers` consists of `n + 1` set functions of the same length `2^n`.
///
/// # Complexity
/// - Time: O(2^n n^2)
/// - Space: O(2^n n)
///
/// # Panics
/// Panics if `layers` is empty.
/// Panics if the length of the set functions is not a power of two.
pub fn ranked_subset_mobius<R: Ring>(ring: &R, mut layers: Vec<Vec<R::Value>>) -> Vec<R::Value> {
    for layer in &mut layers {
        subset_mobius(ring, layer);
    }
    (0..layers[0].len())
        .map(|s| std::mem::replace(&mut layers[s.count_ones() as usize][s], ring.zero()))
        .collect()
}
