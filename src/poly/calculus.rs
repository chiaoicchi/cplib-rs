use crate::algebra::{Field, Semiring};

/// The formal derivative of `f`.
///
/// # Definition
/// The `R`-linear map `x^k -> k x^{k-1}`: the result `g` has `g(k) = (k + 1) f(k + 1)`. It has
/// length `n - 1`.
///
/// # Complexity
/// - Time: O(n)
/// - Space: O(n)
pub fn poly_derivative<R: Semiring>(ring: &R, f: &[R::Value]) -> Vec<R::Value> {
    let one = ring.one();
    let mut k = ring.zero();
    f.iter()
        .skip(1)
        .map(|v| {
            k = ring.add(&k, &one);
            ring.mul(&k, v)
        })
        .collect()
}

/// The integral of `f` vanishing at `0`.
///
/// # Definition
/// The `R`-linear map `x^k -> x^{k+1} / (k + 1)`, the right inverse of the derivative with zero
/// constant term: the result `g` has `g(0) = 0` and `g(k + 1) = f(k) / (k + 1)`. It has length
/// `n + 1`.
///
/// # Contract
/// `1, ..., n` are invertible in `R`,that is the characteristic is `0` or greater than
/// `f.len()`.
///
/// # Complexity
/// - Time: O(n)
/// - Space: O(n)
pub fn poly_integral<R: Field>(ring: &R, f: &[R::Value]) -> Vec<R::Value> {
    if f.is_empty() {
        return Vec::new();
    }
    let one = ring.one();
    let minus_one = ring.neg(&one);
    let mut fact = Vec::with_capacity(f.len() + 1);
    fact.push(ring.one());
    let mut k = ring.zero();
    for i in 1..=f.len() {
        k = ring.add(&k, &one);
        fact.push(ring.mul(&fact[i - 1], &k));
    }
    let mut acc = ring.inv(&fact[f.len()]);
    let mut g: Vec<R::Value> = (0..=f.len()).map(|_| ring.zero()).collect();
    for i in (0..f.len()).rev() {
        g[i + 1] = ring.mul(&f[i], &ring.mul(&fact[i], &acc));
        acc = ring.mul(&acc, &k);
        k = ring.add(&k, &minus_one);
    }
    g
}
