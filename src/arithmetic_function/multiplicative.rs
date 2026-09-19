use crate::algebra::Monoid;
use crate::num::lpf::Lpf;

/// The value of a multiplicative function from a factorization.
///
/// # Definition
/// A function `f` on the nonzero elements of a unique factorization domain, up to units, into a
/// monoid `M` is multiplicative if `f(1) = 1` and `f(ab) = f(a) f(b)` whenever `a` and `b` are
/// coprime. By unique factorization it is determined by its values at prime powers:
/// `f(m) = Π f(p^e)` for `m = Π p^e`, the empty product being `f(1) = 1`.
///
/// # Contract
/// `factors` is a factorization: the `p` are pairwise non-associate primes and every `e >= 1`.
///
/// # Complexity
/// - Time: O(ω(m)) where `ω` is the number of distinct prime factors
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

/// The table of a multiplicative function on `[1, n]`.
///
/// # Definition
/// `table[m] = f(m) = Π f(p^e)` for `m = Π p^e`, `1 <= m <= n`. `table[0]` is `id()`, `0` being
/// outside the domain.
///
/// # Complexity
/// - Time: O(n log log n)
/// - Space: O(n)
pub fn multiplicative_table<M: Monoid>(
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
        let (mut k, mut e) = (i, 0);
        while k % p == 0 {
            k /= p;
            e += 1;
        }
        let v = monoid.op(&table[k], &f(p, e));
        table.push(v);
    }
    table
}
