//! Lagrangian relaxation of an integer constraint.
//!
//! `F` is an objective on a set `S`, `h` is as integer quantity on `S`, and
//! `f(t) = opt {F(s): s∈S, h(s) = t}` is the value function, where `opt` is `max` or `min` below
//! and `f` is defined on a set `X` of consecutive integers. `f` is known only through `solve`,
//! where `solve(λ)` returns `(opt_{s∈S} (F(s) - λh(s)), h(s*))` for some `s*` attaining the
//! optimum.

/// `f(k)` for a concave `f`, where `opt` is `max`.
///
/// # Contract
/// - `f` is concave and `k` is in `X`.
/// - `solve(λ)` is as in the module documentation for every `λ` in `[lo, hi]`, and for some `λ` in
///   `[lo, hi]`, `h(s) = k` for some `s` attaining the maximum.
/// - No computation in `solve(λ)` or of `λk` overflows `i64` for `λ` in `[lo, hi]`.
///
/// # Complexity
/// - Time: O(log (hi - lo + 1)) calls of `solve`
/// - Space: O(1)
///
/// # Panics
/// Panics if `lo > hi`.
pub fn lagrangian_max(
    k: i64,
    mut lo: i64,
    mut hi: i64,
    mut solve: impl FnMut(i64) -> (i64, i64),
) -> i64 {
    assert!(lo <= hi, "lo must no exceed hi: lo={lo}, hi={hi}");
    let mut eval = |lambda: i64| {
        let (value, t) = solve(lambda);
        (value + lambda * k, t)
    };
    let (mut best, t_lo) = eval(lo);
    if t_lo < k {
        return best;
    }
    let (value, t_hi) = eval(hi);
    best = best.min(value);
    if t_hi >= k {
        return best;
    }
    while hi - lo > 1 {
        let mid = lo + (hi - lo) / 2;
        let (value, t) = eval(mid);
        best = best.min(value);
        if t >= k {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    best
}

/// `f(k)` for a convex `f`, where `opt` is `min`.
///
/// # Contract
/// - `f` is convex and `k` is in `X`.
/// - `solve(λ)` is as in the module documentation for every `λ` in `[lo, hi]`, and for some `λ` in
///   `[lo, hi]`, `h(s) = k` for some `s` attaining the minimum.
/// - No computation in `solve(λ)` or of `λk` overflows `i64` for `λ` in `[lo, hi]`.
///
/// # Complexity
/// - Time: O(log (hi - lo + 1)) calls of `solve`
/// - Space: O(1)
///
/// # Panics
/// Panics if `lo > hi`.
pub fn lagrangian_min(k: i64, lo: i64, hi: i64, mut solve: impl FnMut(i64) -> (i64, i64)) -> i64 {
    -lagrangian_max(k, -hi, -lo, |mu| {
        let (value, t) = solve(-mu);
        (-value, t)
    })
}
