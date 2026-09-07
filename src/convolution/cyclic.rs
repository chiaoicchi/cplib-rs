use crate::algebra::RootOfUnity;
use crate::algebra::cyclic::Cyclic;
use crate::convolution::{InverseTransform, Transform, convolve_by_transform};

impl<R: RootOfUnity> Transform<R> for Cyclic {
    /// The discrete Fourier transform `(Ff)(s) = Σ_x w^{sx} f(x)`, where `w` is a primitive
    /// `n`-th root of unity.
    ///
    /// # Definition
    /// The rows of `F` are the characters `x -> w^{sx}` of `Z/nZ`; they exist iff `R` has a
    /// primitive `n`-th root of unity, hence the bound `RootOfUnity`.
    ///
    /// `R[t]/(t^n - 1) = R[t]/(t^{n/2} - 1) x R[t]/(t^{n/2} + 1)` by the Chinese remainder theorem,
    /// and `t -> wu` carries the second factor onto `R[u]/(u^{n/2} - 1)`; iterating `k` times
    /// factors `F` into `k` sparse stages, each with twiddle factors `w^j` (Cooley-Tukey).
    ///
    /// # Complexity
    /// - Time: O(n log n)
    /// - Space: O(log n)
    ///
    /// # Panics
    /// Panics if `f.len() != self.len()`, if `n` is not a power of two, or if `R` has no primitive
    /// `n`-th root of unity.
    fn transform(&self, ring: &R, f: &mut [R::Value]) {
        let w = root(self, ring, f);
        fft(ring, f, w);
    }
}

impl<R: RootOfUnity> InverseTransform<R> for Cyclic {
    /// The inverse discrete Fourier transform `F^{-1} = n^{-1} F'`, where `F'` uses `w^{-1}`.
    ///
    /// # Definition
    /// `Σ_s w^{s(x-y)} = n [x = y]`, so `F' F = n id`. A primitive `n`-th root of unity exists
    /// only if `n != 0` in `R` (otherwise `t^n - 1` has repeated roots), so `n^{-1}` exists
    /// whenever `F` does.
    ///
    /// # Complexity
    /// - Time: O(n log n)
    /// - Space: O(log n)
    ///
    /// # Panics
    /// Panics if `f.len() != self.len()`, if `n` is not a power of two, or if `R` has no primitive
    /// `n`-th root of unity.
    fn inverse_transform(&self, ring: &R, f: &mut [R::Value]) {
        let w = root(self, ring, f);
        fft(ring, f, ring.inv(&w));
        let mut n = ring.one();
        for _ in 0..f.len().trailing_zeros() {
            n = ring.add(&n, &n);
        }
        let inv_n = ring.inv(&n);
        for x in f.iter_mut() {
            *x = ring.mul(x, &inv_n);
        }
    }
}

/// Returns a primitive `f.len()`-th root of unity in `R`.
///
/// # Panics
/// Panics if `f.len() != c.len()`, if the length is not a power of two, or `R` has no primitive
/// root of unity of that order.
fn root<R: RootOfUnity>(c: &Cyclic, ring: &R, f: &[R::Value]) -> R::Value {
    let n = f.len();
    assert_eq!(n, c.order(), "length must be {}: len={n}", c.order());
    assert!(
        n.is_power_of_two(),
        "length must be a power of two: len={n}"
    );
    ring.root_of_unity(n)
        .unwrap_or_else(|| panic!("no primitive {n}-th root of unity"))
}

/// Computes `(Ff)(s) = Σ_x w^{sx} f(x)` in place.
///
/// # Contract
/// `f.len()` is a power of two and `w` is a primitive `f.len()`-th root of unity.
///
/// # Complexity
/// - Time: O(n log n)
/// - Space: O(log n)
fn fft<R: RootOfUnity>(ring: &R, f: &mut [R::Value], w: R::Value) {
    let n = f.len();
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            f.swap(i, j);
        }
    }
    let mut roots = vec![w];
    while roots.len() < n.trailing_zeros() as usize {
        let last = roots.last().unwrap();
        roots.push(ring.mul(last, last));
    }
    let mut len = 2;
    while len <= n {
        let half = len / 2;
        let w_len = roots.pop().unwrap();
        for start in (0..n).step_by(len) {
            let mut wj = ring.one();
            for j in start..start + half {
                let t = ring.mul(&f[j + half], &wj);
                (f[j], f[j + half]) = (ring.add(&f[j], &t), ring.add(&f[j], &ring.neg(&t)));
                wj = ring.mul(&wj, &w_len);
            }
        }
        len <<= 1;
    }
}

/// Returns the product of `f` and `g` in `R[t]`, of length `f.len() + g.len() - 1`.
///
/// # Definition
/// `R[N] = R[t]` embeds into `R[Z/nZ] = R[t]/(t^n - 1)`, injectively on polynomials of degree less
/// than `n`. For `n >= f.len() + g.len() - 1` the product has degree less than `n`, so it is
/// recovered from the cyclic convolution of length `n` without wrap-around.
///
/// # Complexity
/// - Time: O(n log n) with `n` the least power of two at least `f.len() + g.len() - 1`
///
/// # Panics
/// Panics if `R` has no primitive `n`-th root of unity.
pub fn convolve_poly<R: RootOfUnity>(
    ring: &R,
    mut f: Vec<R::Value>,
    mut g: Vec<R::Value>,
) -> Vec<R::Value> {
    if f.is_empty() || g.is_empty() {
        return vec![];
    }
    let len = f.len() + g.len() - 1;
    let n = len.next_power_of_two();
    f.resize_with(n, || ring.zero());
    g.resize_with(n, || ring.zero());
    let mut h = convolve_by_transform(&Cyclic::new(n), ring, f, g);
    h.truncate(len);
    h
}
