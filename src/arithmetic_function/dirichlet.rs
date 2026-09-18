use crate::algebra::Semiring;

/// The Dirichlet convolutio of `f` and `g`.
///
/// # Definition
/// `(fg)(d) = Σ_{xy=d} f(x) g(y)` for `d < n`, the product of the monoid algebra `R[(N, x)]` with
/// the terms of `xy >= n` discarded; `x * 0 = 0`, so `0` absorbs.
///
/// # Complexity
/// - Time: O(n log n)
/// - Space: O(n)
///
/// # Panics
/// Panics if `f` and `g` differ in length.
pub fn dirichlet_convolve<R: Semiring>(
    ring: &R,
    f: Vec<R::Value>,
    g: Vec<R::Value>,
) -> Vec<R::Value> {
    assert_eq!(
        f.len(),
        g.len(),
        "f and g must have the same length: f={}, g={}",
        f.len(),
        g.len()
    );
    let n = f.len();
    let mut h: Vec<R::Value> = (0..n).map(|_| ring.zero()).collect();
    for a in 1..n {
        for b in 1..=(n - 1) / a {
            h[a * b] = ring.add(&h[a * b], &ring.mul(&f[a], &g[b]));
        }
    }
    if n > 0 {
        for b in 0..n {
            h[0] = ring.add(&h[0], &ring.mul(&f[0], &g[b]));
        }
        for a in 1..n {
            h[0] = ring.add(&h[0], &ring.mul(&f[a], &g[0]));
        }
    }
    h
}
