use crate::algebra::{Field, RootOfUnity};
use crate::cyclic::cyclic_convolve;
use crate::poly::subproduct_tree::SubproductTree;

/// The values at the points `xs` of the polynomial with coefficients `f`.
///
/// # Definition
/// For `n = xs.len()` and `m = f.len()`, with `x(j)` the `j`-th entry of `xs`. Returns
/// `f(x(j)) = Σ_k f(k) x(j)^k` for `j` in `[0, n)`.
///
/// # Complexity
/// - Time: O(m log m + n log^2 n)
/// - Space: O(m + n log n)
///
/// # Panics
/// Panics if `R` has no primitive `N`-th root of unity, where `N` is the least power of two at
/// least `max(n, a)`, with `a = 0` if `n = 0`, `a = m` if `m <= 32`, and `a = 2m - 1` otherwise.
pub fn multipoint_evaluate<R: RootOfUnity + Clone>(
    ring: &R,
    f: &[R::Value],
    xs: &[R::Value],
) -> Vec<R::Value> {
    SubproductTree::new(ring.clone(), xs).evaluate(f)
}

/// The value at `c` of the polynomial of degree less than `n` taking the values `ys` on `[0, n)`.
///
/// # Definition
/// For `n = ys.len()`, with `y(i)` the `i`-th entry of `ys`, let `f` in `R[x]` be the unique
/// polynomial of degree less than `n` with `f(i) = y(i)` for `i` in `[0, n)`. Returns
/// `f(c) = Σ_i y(i) Π_{j!=i} (c - j) / Π_{j!=i} (i - j)`. For `n = 0`, `f` is zero.
///
/// # Contract
/// `1, ..., n - 1` are invertible in `R`, that is the characteristic is `0` or greater than
/// `n - 1`. Under it `f` exists and is unique, as the Vandermonde determinant `Π_{k<n} k!` is
/// invertible.
///
/// # Complexity
/// - Time: O(n)
/// - Space: O(n)
pub fn iota_evaluate<R: Field>(ring: &R, ys: &[R::Value], c: &R::Value) -> R::Value {
    let n = ys.len();
    if n == 0 {
        return ring.zero();
    }
    let one = ring.one();
    let minus_one = ring.neg(&one);

    let mut k = ring.zero();
    let mut fact = ring.one();
    for _ in 1..n {
        k = ring.add(&k, &one);
        fact = ring.mul(&fact, &k);
    }

    let mut inv_fact: Vec<R::Value> = (0..n).map(|_| ring.one()).collect();
    let mut suf: Vec<R::Value> = (0..=n).map(|_| ring.one()).collect();
    let mut inv = ring.inv(&fact);
    for i in (0..n).rev() {
        suf[i] = ring.mul(&suf[i + 1], &ring.add(c, &ring.neg(&k)));
        let next = ring.mul(&inv, &k);
        inv_fact[i] = inv;
        inv = next;
        k = ring.add(&k, &minus_one);
    }

    let mut k = ring.zero();
    let mut pre = ring.one();
    let mut acc = [ring.zero(), ring.zero()];
    for i in 0..n {
        let numer = ring.mul(&pre, &suf[i + 1]);
        let denom = ring.mul(&inv_fact[i], &inv_fact[n - 1 - i]);
        let term = ring.mul(&ring.mul(&ys[i], &numer), &denom);
        let s = (n - 1 - i) & 1;
        acc[s] = ring.add(&acc[s], &term);
        pre = ring.mul(&pre, &ring.add(c, &ring.neg(&k)));
        k = ring.add(&k, &one);
    }
    ring.add(&acc[0], &ring.neg(&acc[1]))
}

/// The values `f(c), f(c + 1), ..., f(c + m - 1)` of the polynomial of degree less than `n` taking
/// the values `ys` on `[0, n)`.
///
/// # Definition
/// For `n = ys.len()`, with `y(i)` the `i`-th entry of `ys`, let `f` in `R[x]` be the unique
/// polynomial of degree less than `n` with `f(i) = y(i)` for `i` in `[0, n)`. Returns `f(c + k)`
/// for `k` in `[0, m)`. For `m = 1` it is [`iota_evaluate`].
///
/// # Contract
/// `1, ..., n + m - 2` are invertible in `R`.
///
/// # Complexity
/// - Time: O((n + m) log (n + m))
/// - Space: O(n + m)
///
/// # Panics
/// Panics if `n > 0`, `m > 0`, and `R` has no primitive `N`-th root of unity, where `N` is the
/// least power of two at least `n + m - 1`.
pub fn iota_shift<R: RootOfUnity<Value: PartialEq>>(
    ring: &R,
    ys: &[R::Value],
    c: &R::Value,
    m: usize,
) -> Vec<R::Value> {
    let n = ys.len();
    if m == 0 {
        return Vec::new();
    }
    if n == 0 {
        return (0..m).map(|_| ring.zero()).collect();
    }
    let zero = ring.zero();
    let one = ring.one();
    let minus_one = ring.neg(&one);
    let len = n + m - 1;
    let size = len.next_power_of_two();

    let mut k = ring.zero();
    let mut fact = ring.one();
    for _ in 1..n {
        k = ring.add(&k, &one);
        fact = ring.mul(&fact, &k);
    }

    let mut x: Vec<R::Value> = Vec::with_capacity(n + m - 1);
    let mut u0 = None;
    let mut xu = ring.add(c, &ring.neg(&k));
    for u in 0..len {
        let next = ring.add(&xu, &one);
        if xu == zero {
            u0 = Some(u);
            x.push(ring.one());
        } else {
            x.push(xu);
        }
        xu = next;
    }

    let mut inv_fact: Vec<R::Value> = (0..n).map(|_| ring.one()).collect();
    let mut inv = ring.inv(&fact);
    for i in (0..n).rev() {
        let next = ring.mul(&inv, &k);
        inv_fact[i] = inv;
        inv = next;
        k = ring.add(&k, &minus_one);
    }
    let mut a: Vec<R::Value> = Vec::with_capacity(size);
    for i in 0..n {
        let t = ring.mul(&ys[i], &ring.mul(&inv_fact[i], &inv_fact[n - 1 - i]));
        a.push(if (n - 1 - i) & 1 == 1 {
            ring.neg(&t)
        } else {
            t
        });
    }

    let mut q: Vec<R::Value> = Vec::with_capacity(len);
    q.push(ring.mul(&x[0], &one));
    for j in 1..len {
        q.push(ring.mul(&q[j - 1], &x[j]));
    }
    let mut inv_q: Vec<R::Value> = (0..len).map(|_| ring.one()).collect();
    let mut inv = ring.inv(&q[len - 1]);
    for j in (0..len).rev() {
        let next = ring.mul(&inv, &x[j]);
        inv_q[j] = inv;
        inv = next;
    }

    let mut b: Vec<R::Value> = Vec::with_capacity(size);
    for u in 0..len {
        b.push(if u0 == Some(u) {
            ring.zero()
        } else if u == 0 {
            ring.mul(&inv_q[0], &one)
        } else {
            ring.mul(&q[u - 1], &inv_q[u])
        });
    }

    a.resize_with(size, || ring.zero());
    b.resize_with(size, || ring.zero());
    let s = cyclic_convolve(ring, a, b);

    let mut values: Vec<R::Value> = Vec::with_capacity(m);
    for k in 0..m {
        let window = if k == 0 {
            ring.mul(&q[n - 1], &one)
        } else {
            ring.mul(&q[k + n - 1], &inv_q[k - 1])
        };
        values.push(ring.mul(&window, &s[k + n - 1]));
    }
    if let Some(u0) = u0 {
        for (k, value) in values.iter_mut().enumerate() {
            if let Some(i) = (k + n - 1).checked_sub(u0).filter(|&i| i < n) {
                *value = ring.add(&ys[i], &zero);
            }
        }
    }
    values
}
