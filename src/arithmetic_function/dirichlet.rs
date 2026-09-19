use crate::algebra::Semiring;

/// The Dirichlet convolution of `f` and `g`.
///
/// # Definition
/// `(fg)(d) = Σ_{xy=d} f(x) g(y)` for `1 <= d < n`, and `(fg)(0)` collects the terms with
/// `xy = 0` or `xy >= n`: the product of the monoid algebra `R[(N, x)/I]`, where `I = {0} ∪ [n, ∞)`
/// is an ideal of `(N, x)` and `0` its class.
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
        let fa = &f[a];
        for (hd, gb) in h[a..].iter_mut().step_by(a).zip(&g[1..]) {
            *hd = ring.add(hd, &ring.mul(fa, gb));
        }
    }
    let mut suffix: Vec<R::Value> = (0..=n).map(|_| ring.zero()).collect();
    for b in (0..n).rev() {
        suffix[b] = ring.add(&suffix[b + 1], &g[b]);
    }
    for a in 1..n {
        let overflow = ring.add(&g[0], &suffix[n.div_ceil(a)]);
        h[0] = ring.add(&h[0], &ring.mul(&f[a], &overflow));
    }
    if n > 0 {
        h[0] = ring.add(&h[0], &ring.mul(&f[0], &suffix[0]));
    }
    h
}
