use crate::algebra::RootOfUnity;
use crate::cyclic::{dif, dit, twiddles};
use crate::fps::fps_convolve;
use crate::poly::calculus::{poly_derivative, poly_integral};

/// The inverse of `f` in `R[[x]]/(x^n)`, as its `n` coefficients.
///
/// # Definition
/// `g(0) = f(0)^{-1}` and `g(k) = -f(0)^{-1} Σ_{i=1,...,k} f(i) g(k - i)` for `k >= 1`.
///
/// # Contract
/// `f(0)` is invertible in `R`, unless `n = 0`.
///
/// # Complexity
/// - Time: O(n log n)
/// - Space: O(n)
///
/// # Panics
/// Panics if `R` has no primitive `N`-th root of unity, where `N` is the least power of two at
/// least `n`.
pub fn fps_inv<R: RootOfUnity>(ring: &R, f: &[R::Value], n: usize) -> Vec<R::Value> {
    if n == 0 {
        return Vec::new();
    }
    let zero = ring.zero();
    let mut g = vec![ring.inv(&f[0])];
    if n == 1 {
        return g;
    }
    let size = n.next_power_of_two();
    let omega = ring
        .root_of_unity(size)
        .unwrap_or_else(|| panic!("no primitive {size}-th root of unity"));
    let inv_omega = ring.inv(&omega);
    let table = twiddles(ring, omega, size);
    let inv_table = twiddles(ring, inv_omega, size);
    let half = ring.inv(&ring.add(&ring.one(), &ring.one()));
    let mut inv_len = ring.one();
    let mut m = 1;
    while m < n {
        let len = m << 1;
        inv_len = ring.mul(&inv_len, &half);
        let mut e: Vec<R::Value> = f.iter().take(len).map(|x| ring.add(x, &zero)).collect();
        e.resize_with(len, || ring.zero());
        let mut gd: Vec<R::Value> = g.iter().map(|x| ring.add(x, &zero)).collect();
        gd.resize_with(len, || ring.zero());
        dif(ring, &mut gd, &table);
        dif(ring, &mut e, &table);
        for (x, y) in e.iter_mut().zip(&gd) {
            *x = ring.mul(x, y);
        }
        dit(ring, &mut e, &inv_table);
        for x in e.iter_mut().take(m) {
            *x = ring.zero();
        }
        dif(ring, &mut e, &table);
        for (x, y) in e.iter_mut().zip(&gd) {
            *x = ring.mul(x, y);
        }
        dit(ring, &mut e, &inv_table);
        let scale = ring.neg(&ring.mul(&inv_len, &inv_len));
        g.extend(e[m..].iter().map(|x| ring.mul(x, &scale)));
        m = len;
    }
    g.truncate(n);
    g
}

/// The logarithm `g` of `f` in `R[[x]]/(x^n)`, with `g(0) = 0` and `g' = f' / f`, as its `n`
/// coefficients.
///
/// # Definition
/// `log f = Σ_{k>=1} (-1)^{k+1} (f - 1)^k / k`.
///
/// # Contract
/// `f(0) = 1`, unless `n = 0`. `1, ..., n - 1` are invertible in `R`.
///
/// # Complexity
/// - Time: O(n log n)
/// - Space: O(n)
///
/// # Panics
/// Panics only if `R` has no primitive `N`-th root of unity, where `N` is the least power of two at
/// least `2n`.
pub fn fps_log<R: RootOfUnity>(ring: &R, f: &[R::Value], n: usize) -> Vec<R::Value> {
    if n == 0 {
        return Vec::new();
    }
    let h = fps_convolve(
        ring,
        poly_derivative(ring, &f[..f.len().min(n)]),
        fps_inv(ring, f, n - 1),
        n - 1,
    );
    poly_integral(ring, &h)
}

/// The exponential `g` of `f` in `R[[x]]/(x^n)`, with `g(0) = 1` and `g' = f'g`, as its `n`
/// coefficients.
///
/// # Definition
/// `exp f = Σ_{k>=0} f^k / k!`.
///
/// # Contract
/// `f(0) = 0`, unless `n = 0`. `1, ..., n - 1` are invertible in `R`.
///
/// # Complexity
/// - Time: O(n log n)
/// - Space: O(n)
///
/// # Panics
/// Panics only if `R` has no primitive `N`-th root of unity, where `N` is the least power of two at
/// least `2n`.
pub fn fps_exp<R: RootOfUnity>(ring: &R, f: &[R::Value], n: usize) -> Vec<R::Value> {
    if n == 0 {
        return Vec::new();
    }
    let mut g = vec![ring.one()];
    let mut m = 1;
    while m < n {
        m = (m << 1).min(n);
        let mut t: Vec<R::Value> = fps_log(ring, &g, m).iter().map(|x| ring.neg(x)).collect();
        t[0] = ring.add(&t[0], &ring.one());
        for (t, x) in t.iter_mut().zip(f) {
            *t = ring.add(t, x);
        }
        g = fps_convolve(ring, g, t, m);
    }
    g
}

/// The power `f^k` in `R[[x]]/(x^n)`, as its `n` coefficients.
///
/// # Definition
/// `f^0 = 1`, also for `f = 0`, and `f^k = f f^{k-1}` for `k >= 1`.
///
/// # Contract
/// If `k >= 1` and `kv < n`, `1, ..., n - v - 1` are invertible in `R`, where `v` is the least
/// index with `f(v) != 0`, if any.
///
/// # Complexity
/// - Time: O(n log n + log k)
/// - Space: O(n)
///
/// # Panics
/// Panics only if `R` has no primitive `N`-th root of unity, where `N` is the least power of two at
/// least `2n`.
pub fn fps_pow<R: RootOfUnity<Value: PartialEq>>(
    ring: &R,
    f: &[R::Value],
    n: usize,
    k: u64,
) -> Vec<R::Value> {
    let zero = ring.zero();
    let one = ring.one();
    let mut g: Vec<R::Value> = (0..n).map(|_| ring.zero()).collect();
    if n == 0 {
        return g;
    }
    if k == 0 {
        g[0] = one;
        return g;
    }
    let Some(v) = f.iter().take(n).position(|c| *c != zero) else {
        return g;
    };
    if v > 0 && k >= (n as u64).div_ceil(v as u64) {
        return g;
    }

    let c = &f[v];
    let mut k_image = ring.zero();
    let mut c_k = ring.one();
    for i in (0..u64::BITS - k.leading_zeros()).rev() {
        k_image = ring.add(&k_image, &k_image);
        c_k = ring.mul(&c_k, &c_k);
        if k >> i & 1 == 1 {
            k_image = ring.add(&k_image, &one);
            c_k = ring.mul(&c_k, c);
        }
    }

    let c_inv = ring.inv(c);
    let u: Vec<R::Value> = f[v..]
        .iter()
        .take(n - v)
        .map(|x| ring.mul(x, &c_inv))
        .collect();
    let l: Vec<R::Value> = fps_log(ring, &u, n - v)
        .iter()
        .map(|x| ring.mul(x, &k_image))
        .collect();
    let w = fps_exp(ring, &l, n - v);
    let shift = k as usize * v;
    for (g, w) in g[shift..].iter_mut().zip(&w) {
        *g = ring.mul(w, &c_k);
    }
    g
}
