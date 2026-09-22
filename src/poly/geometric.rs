use crate::algebra::Field;
use crate::poly::evaluation::iota_evaluate;

/// The value of the series `Σ_{i>=0} r^i f(i)`, for a polynomial `f` given by its values.
///
/// # Definition
/// For `m = ys.len()`, let `f` in `R[x]` be the polynomial of degree less than `m` with
/// `f(i) = y(i)` for `i` in `[0, m)`. The generating series `A(x) = Σ_{i>=0} r^i f(i) x^i` is the
/// rational function `P(rx) / (1 - rx)^m` with `deg P < m`. Returns `A(1)`, which is the sum of the
/// series whenever that converges.
///
/// Equivalently, `A(1)` is the unique constant `C` such that `Σ_{i<k} r^i f(i) = C + r^k h(k)` for
/// all `k >= 0`, for some polynomial `h` of degree less than `m`.
///
/// # Contract
/// `1, ..., m - 1` are invertible in `R`, `1 - r` is invertible in `R`, unless `m = 0`.
///
/// # Complexity
/// - Time: O(m)
/// - Space: O(m)
pub fn iota_geometric_series<R: Field>(ring: &R, ys: &[R::Value], r: &R::Value) -> R::Value {
    let m = ys.len();
    if m == 0 {
        return ring.zero();
    }
    let one = ring.one();
    let minus_one = ring.neg(&one);

    let mut k = ring.zero();
    let mut fact = ring.one();
    for _ in 1..=m {
        k = ring.add(&k, &one);
        fact = ring.mul(&fact, &k);
    }
    let mut inv_fact: Vec<R::Value> = (0..m).map(|_| ring.one()).collect();
    let mut inv = ring.inv(&fact);
    for i in (0..m).rev() {
        let next = ring.mul(&inv, &k);
        inv_fact[i] = inv;
        inv = next;
        k = ring.add(&k, &minus_one);
    }

    let mut pow_r: Vec<R::Value> = Vec::with_capacity(m + 1);
    pow_r.push(ring.one());
    for i in 1..=m {
        pow_r.push(ring.mul(&pow_r[i - 1], r));
    }
    let mut s = Vec::with_capacity(m + 1);
    let mut g = ring.zero();
    for k in 0..=m {
        let t = ring.mul(&pow_r[m - k], &g);
        s.push(if (m - k) & 1 == 1 { ring.neg(&t) } else { t });
        if k < m {
            g = ring.add(&g, &ring.mul(&pow_r[k], &ys[k]));
        }
    }

    let mut sum = ring.zero();
    for j in 0..m {
        let binom = ring.mul(&fact, &ring.mul(&inv_fact[j], &inv_fact[m - 1 - j]));
        sum = ring.add(&sum, &ring.mul(&binom, &ring.add(&s[j], &s[j + 1])));
    }
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
/// For `f` as in [`iota_geometric_series`], with `r^0 = 1` also for `r = 0`. For `f = 1` it is
/// [`geometric_sum`](crate::algebra::power::geometric_sum).
///
/// # Contract
/// `1, ..., m - 1` are invertible in `R`.
///
/// # Complexity
/// - Time: O(m + log n)
/// - Space: O(m)
pub fn iota_geometric_sum<R: Field<Value: PartialEq>>(
    ring: &R,
    ys: &[R::Value],
    r: &R::Value,
    n: u64,
) -> R::Value {
    let m = ys.len();
    if m == 0 || n == 0 {
        return ring.zero();
    }
    let zero = ring.zero();
    let one = ring.one();
    if *r == zero {
        return ring.add(&ys[0], &zero);
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
        let mut g = Vec::with_capacity(m + 1);
        g.push(ring.zero());
        for i in 0..m {
            g.push(ring.add(&g[i], &ys[i]));
        }
        let mut m_image = ring.zero();
        for i in (0..u64::BITS - m.leading_zeros()).rev() {
            m_image = ring.add(&m_image, &m_image);
            if m >> i & 1 == 1 {
                m_image = ring.add(&m_image, &one);
            }
        }
        if m_image == zero {
            let q = n / m as u64;
            let mut q_image = ring.zero();
            for i in (0..u64::BITS - m.leading_zeros()).rev() {
                q_image = ring.add(&q_image, &q_image);
                if q >> i & 1 == 1 {
                    q_image = ring.add(&q_image, &one);
                }
            }
            return ring.add(&ring.mul(&q_image, &g[m]), &g[(n % m as u64) as usize]);
        }
        return iota_evaluate(ring, &g, &n_image);
    }

    let c = iota_geometric_series(ring, ys, r);
    let minus_c = ring.neg(&c);
    let inv_r = ring.inv(r);
    let mut h = Vec::with_capacity(m);
    let mut g = ring.zero();
    let mut pow_r = ring.one();
    let mut pow_inv_r = ring.one();
    for y in ys {
        h.push(ring.mul(&ring.add(&g, &minus_c), &pow_inv_r));
        g = ring.add(&g, &ring.mul(&pow_r, y));
        pow_r = ring.mul(&pow_r, r);
        pow_inv_r = ring.mul(&pow_inv_r, &inv_r);
    }
    let h_n = iota_evaluate(ring, &h, &n_image);
    ring.add(&c, &ring.mul(&r_n, &h_n))
}
