//! Multiplicative functions.
//!
//! A function `g` on the positive integers into a monoid `M` is multiplicative if `g(1) = id()`
//! and `g(ab) = g(a) g(b)` whenever `a` and `b` are coprime. Throughout this module `g` is the
//! multiplicative function with `g(p^e) = f(p, e)` for every prime power `p^e`.

use crate::algebra::{Commutative, Monoid};
use crate::arithmetic::lpf::Lpf;

/// The value of `g` from a factorization.
///
/// # Definition
/// `Π f(p, e)` over the pairs `(p, e)` of `factors`, in this order, and `id()` if `factors` is
/// empty. For the factorization `m = Π p^e`, it is `g(m)`.
///
/// # Complexity
/// - Time: O(k), where `k = factors.len()`
/// - Space: O(1)
pub fn multiplicative<M: Monoid, P>(
    monoid: &M,
    factors: &[(P, u32)],
    f: impl Fn(&P, u32) -> M::Value,
) -> M::Value {
    factors
        .iter()
        .fold(monoid.id(), |acc, (p, e)| monoid.op(&acc, &f(p, *e)))
}

/// The table of `g` on `[1, n]`.
///
/// # Definition
/// `table[m] = g(m)` for `1 <= m <= n`, so `table[1] = id()`. `table[0] = id()`.
///
/// # Complexity
/// - Time: O(n log log n)
/// - Space: O(n)
pub fn multiplicative_table<M: Monoid + Commutative>(
    monoid: &M,
    n: usize,
    f: impl Fn(usize, u32) -> M::Value,
) -> Vec<M::Value> {
    let lpf = Lpf::new(n);
    let mut table = Vec::with_capacity(n + 1);
    for i in 0..=n {
        if i < 2 {
            table.push(monoid.id());
            continue;
        }
        let p = lpf.lpf(i);
        let (mut k, mut e, mut q) = (i / p, 0, p);
        while k % p == 0 {
            k /= p;
            e += 1;
            q *= p;
        }
        let v = if k == 1 {
            f(p, e)
        } else {
            monoid.op(&table[k], &table[q])
        };
        table.push(v);
    }
    table
}
