use crate::algebra::RootOfUnity;
use crate::poly::poly_convolve;

/// The polynomial `f(x + c)`.
///
/// # Definition
/// The image of `f` under the `R`-algebra automorphism of `R[x]` sending `x` to `x + c`:
/// `[x^k] f(x + c) = Σ_{i>=k} binom(i, k) f(i) c^{i-k}`. The result has length `n`.
///
/// # Contract
/// `1, ..., n - 1` are invertible in `R`, that is the characteristic is `0` or greater than
/// `n - 1`.
///
/// # Complexity
/// - Time: O(n log n)
/// - Space: O(n)
///
/// # Panics
/// Panics if `n > 32` and `R` has no primitive `N`-th root of unity, where `N` is the least power
/// of two at least `2n - 1`.
pub fn poly_taylor_shift<R: RootOfUnity>(ring: &R, f: &[R::Value], c: &R::Value) -> Vec<R::Value> {
    let n = f.len();
    if n == 0 {
        return Vec::new();
    }
    let one = ring.one();
    let minus_one = ring.neg(&one);

    let mut g: Vec<R::Value> = Vec::with_capacity(n);
    let mut fact = ring.one();
    let mut k = ring.zero();
    for (i, f) in f.iter().enumerate() {
        if i > 0 {
            k = ring.add(&k, &one);
            fact = ring.mul(&fact, &k);
        }
        g.push(ring.mul(&fact, f));
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

    let mut p = poly_convolve(ring, g, h);
    p.truncate(n);
    for (p, inv) in p.iter_mut().zip(inv_fact.iter().rev()) {
        *p = ring.mul(p, inv);
    }
    p.reverse();
    p
}
