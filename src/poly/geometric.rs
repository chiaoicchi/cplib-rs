use crate::algebra::Field;
use crate::poly::interpolation::lagrange_interpolate_iota;

/// The value of the series `Σ_{i>=0} r^i f(i)`, for a polynomial `f` given by its values.
///
/// # Definition
/// For `m = y.len()`, let `f` in `R[x]` be the polynomial of degree less than `m` with
/// `f(i) = y(i)` for `i` in `[0, m)`. The generating series `A(x) = Σ_{i>=0} r^i f(i) x^i` is the
/// rational function `P(rx) / (1 - rx)^m` with `deg P < m`. Returns `A(1)`, which is the sum of the
/// series whenever that converges.
///
/// Equivalently, there are a unique constant `C` and a unique polynomial `h` of degree less than `m`
/// with `Σ_{i<k} r^i f(i) = C + r^k h(k)` for all `k >= 0`, and `A(1) = C`.
///
/// # Contract
/// `1 - r` and `m!` are invertible in `R`.
///
/// # Complexity
/// - Time: O(m)
/// - Space: O(m)
pub fn poly_geometric_series<R: Field>(ring: &R, y: &[R::Value], r: &R::Value) -> R::Value {
    let m = y.len();
    if m == 0 {
        return ring.zero();
    }
    let one = ring.one();
    let minus_one = ring.neg(&one);

    let mut k = ring.zero();
    let mut fact = ring.one();
    let mut pow_r: Vec<R::Value> = Vec::with_capacity(m + 1);
    pow_r.push(ring.one());
    for i in 1..=m {
        k = ring.add(&k, &one);
        fact = ring.mul(&fact, &k);
        pow_r.push(ring.mul(&pow_r[i - 1], r));
    }

    let mut inv_fact: Vec<R::Value> = (0..=m).map(|_| ring.one()).collect();
    let mut inv = ring.inv(&fact);
    for i in (0..=m).rev() {
        let next = ring.mul(&inv, &k);
        inv_fact[i] = inv;
        inv = next;
        k = ring.add(&k, &minus_one);
    }

    let mut g = ring.zero();
    let mut acc = [ring.zero(), ring.zero()];
    for k in 0..=m {
        let binom = ring.mul(&inv_fact[k], &inv_fact[m - k]);
        let term = ring.mul(&binom, &ring.mul(&pow_r[m - k], &g));
        let s = (m - k) & 1;
        acc[s] = ring.add(&acc[s], &term);
        if k < m {
            g = ring.add(&g, &ring.mul(&pow_r[k], &y[k]));
        }
    }
    let sum = ring.mul(&fact, &ring.add(&acc[0], &ring.neg(&acc[1])));
    let base = ring.add(&one, &ring.neg(r));
    let mut denom = ring.one();
    for i in (0..usize::BITS - m.leading_zeros()).rev() {
        denom = ring.mul(&denom, &denom);
        if m >> i & 1 == 1 {
            denom = ring.mul(&denom, &base);
        }
    }
    ring.mul(&sum, &ring.inv(&denom))
}

/// The sum `Σ_{i=0,...,n-1} r^i f(i)`, for a polynomial `f` given by its values.
///
/// # Definition
/// For `f` as in [`poly_geometric_series`], with `r^0 = 1` also for `r = 0`. For `f = 1` it is
/// [`geometric_sum`](crate::algebra::power::geometric_sum).
///
/// # Contract
/// `m!` is invertible in `R`.
///
/// # Complexity
/// - Time: O(m + log n)
/// - Space: O(m)
pub fn poly_geometric_sum<R: Field<Value: PartialEq>>(
    ring: &R,
    y: &[R::Value],
    r: &R::Value,
    n: u64,
) -> R::Value {
    let m = y.len();
    if m == 0 || n == 0 {
        return ring.zero();
    }
    let zero = ring.zero();
    let one = ring.one();
    if *r == zero {
        return ring.add(&y[0], &zero);
    }

    let mut n_image = ring.zero();
    let mut r_n = ring.one();
    for i in (0..u64::BITS - n.leading_zeros()).rev() {
        n_image = ring.add(&n_image, &n_image);
        r_n = ring.mul(&r_n, &r_n);
        if n >> i & 1 == 1 {
            n_image = ring.add(&n_image, &one);
            r_n = ring.mul(&r_n, r);
        }
    }

    if *r == one {
        let mut g: Vec<R::Value> = Vec::with_capacity(m + 1);
        g.push(ring.zero());
        for i in 0..m {
            g.push(ring.add(&g[i], &y[i]));
        }
        return lagrange_interpolate_iota(ring, &g, &n_image);
    }

    let c = poly_geometric_series(ring, y, r);
    let minus_c = ring.neg(&c);
    let inv_r = ring.inv(r);
    let mut h: Vec<R::Value> = Vec::with_capacity(m);
    let mut g = ring.zero();
    let mut pow_r = ring.one();
    let mut pow_inv_r = ring.one();
    for yk in y {
        h.push(ring.mul(&ring.add(&g, &minus_c), &pow_inv_r));
        g = ring.add(&g, &ring.mul(&pow_r, yk));
        pow_r = ring.mul(&pow_r, r);
        pow_inv_r = ring.mul(&pow_inv_r, &inv_r);
    }
    let h_n = lagrange_interpolate_iota(ring, &h, &n_image);
    ring.add(&c, &ring.mul(&r_n, &h_n))
}
