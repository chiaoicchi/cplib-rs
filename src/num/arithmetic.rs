use crate::algebra::canonical::Canonical;
use crate::algebra::multiplicative::Multiplicative;
use crate::algebra::{Monoid, Semiring};
use crate::num::prime::spf;

/// Evaluates a multiplicative function from its values `f(p, e) = f(p^e)` at prime powers.
///
/// # Definition
/// A function `f` on the positive integers into a monoid `M` is multiplicative if
/// `f(1) = 1` and `f(ab) = f(a)f(b)` whenever `gcd(a, b) = 1`. By unique factorization such `f` is
/// determined by its values at prime powers: `f(n) = Π f(p^e)` for `n = Π p^e`, the empty product
/// being `f(1) = 1`.
///
/// # Contract
/// `factors` is a factorization: the `p` are distinct primes and every `e >= 1`.
///
/// # Complexity
/// - Time: O(ω(n)) where `ω(n)` is the number of distinct prime factors
/// - Space: O(1)
pub fn multiplicative<M: Monoid>(
    monoid: &M,
    factors: &[(u64, u32)],
    f: impl Fn(u64, u32) -> M::Value,
) -> M::Value {
    factors
        .iter()
        .fold(monoid.id(), |acc, &(p, e)| monoid.op(&acc, &f(p, e)))
}

/// Tabulates a multiplicative function on `[1, n]` from its values `f(p, e) = f(p^e)` at prime
/// powers.
///
/// # Definition
/// `table[m] = f(m) = Π f(p^e)` over `m = Π p^e` for `1 <= m <= n`; in particular
/// `table[1] = 1`. `table[0]` is set to `id()`, `0` being outside the domain.
///
/// # Complexity
/// - Time: O(n log log n)
/// - Space: O(n)
pub fn multiplicative_table<M: Monoid>(
    monoid: &M,
    n: usize,
    f: impl Fn(usize, u32) -> M::Value,
) -> Vec<M::Value>
where
    M::Value: Clone,
{
    let spf = spf(n);
    let mut table = Vec::with_capacity(n + 1);
    for (m, spf) in spf.iter().enumerate() {
        if m < 2 {
            table.push(monoid.id());
            continue;
        }
        let (mut k, mut e) = (m, 0);
        while k % spf == 0 {
            k /= spf;
            e += 1;
        }
        let v = monoid.op(&table[k], &f(*spf, e));
        table.push(v);
    }
    table
}

/// The Dirichlet convolution of two arithmetic functions tabulated on `[1, n]`.
///
/// # Definition
/// `(f * g)(m) = Σ_{ab=m} f(a)g(b)` for `1 <= m <= n`, the sum over ordered pairs of positive
/// integers. Only values of `f` and `g` on `[1, n]` enter, since `ab = m <= n` forces `a, b <= n`,
/// so `{1, ..., n}` is closed under `*`.
/// Index `0` is outside the domain and is set to `zero()`.
///
/// # Complexity
/// - Time: O(n log n)
/// - Space: O(n)
///
/// # Panics
/// Panics if `f.len() != g.len()`.
/// Panics if `f` or `g` is empty.
pub fn dirichlet_convolve<R: Semiring>(ring: &R, f: &[R::Value], g: &[R::Value]) -> Vec<R::Value> {
    assert_eq!(
        f.len(),
        g.len(),
        "length must agree: {} != {}",
        f.len(),
        g.len()
    );
    assert!(f.is_empty(), "f must not be empty");
    let n = f.len() - 1;
    let mut h: Vec<R::Value> = (0..=n).map(|_| ring.zero()).collect();
    for a in 1..=n {
        for b in 1..=n / a {
            h[a * b] = ring.add(&h[a * b], &ring.mul(&f[a], &g[b]));
        }
    }
    h
}

/// Euler's totient `φ(n)` from a factorization of `n`.
///
/// # Definition
/// `φ(n) = #{s in [1, n]: gcd(s, n) = 1} = |(Z/nZ)^*|`. It is multiplicative, since the Chinese
/// remainder theorem gives `(Z/abZ)^* = (Z/aZ)^* x (Z/bZ)^*` for coprime `a`, `b`, with
/// `φ(p^e) = p^e - p^{e-1}`, the multiples of `p` being the residues not coprime to `p^e`. Hence
/// `φ(n) = Π p^{e-1}(p - 1) = n Π (1 - 1/p)`.
///
/// # Contract
/// `factors` is a factorization.
///
/// # Complexity
/// - Time: O(ω(n) log n)
/// - Space: O(1)
pub fn euler_phi(factors: &[(u64, u32)]) -> u64 {
    multiplicative(&Multiplicative(Canonical::new()), factors, |p, e| {
        p.pow(e - 1) * (p - 1)
    })
}

/// Euler's totient `φ(m)` for `1 <= m <= n`.
///
/// # Definition
/// See [`euler_phi`]; tabulated by [`multiplicative_table`].
///
/// # Complexity
/// - Time: O(n log log n)
/// - Space: O(n)
pub fn phi_table(n: usize) -> Vec<u64> {
    multiplicative_table(&Multiplicative(Canonical::new()), n, |p, e| {
        (p as u64).pow(e - 1) * (p as u64 - 1)
    })
}

/// The Möbius function `μ(n)` from a factorization of `n`.
///
/// # Definition
/// `μ(n) = (-1)^r` if `n` is a product of `r` distinct primes and `μ(n) = 0` if `n` has a square
/// factor; in particular `μ(1) = 1`. It is multiplicative with `μ(p) = -1`
/// and `μ(p^e) = 0` for `e >= 2`.
///
/// # Contract
/// `factors` is a factorization.
///
/// # Complexity
/// - Time: O(ω(n))
/// - Space: O(1)
pub fn mobius(factors: &[(u64, u32)]) -> i64 {
    multiplicative(&Multiplicative(Canonical::new()), factors, |_, e| {
        if e == 1 { -1 } else { 0 }
    })
}

/// The Möbius function `μ(m)` for `1 <= m <= n`.
///
/// # Definition
/// See [`mobius`]; tabulated by [`multiplicative_table`].
///
/// # Complexity
/// - Time: O(n log log n)
/// - Space: O(n)
pub fn mobius_table(n: usize) -> Vec<i64> {
    multiplicative_table(&Multiplicative(Canonical::new()), n, |_, e| {
        if e == 1 { -1 } else { 0 }
    })
}

/// The number of divisors `d(n)` from a factorization of `n`.
///
/// # Definition
/// `d(n) = #{m: m|n}`. A divisor of `n = Π p^e` is `Π p^f` with `0 <= f <= e`, so `d` is
/// multiplicative with `d(p^e) = e + 1`.
///
/// # Contract
/// `factors` is a factorization.
///
/// # Complexity
/// - Time: O(ω(n))
/// - Space: O(1)
pub fn num_divisors(factors: &[(u64, u32)]) -> u64 {
    multiplicative(&Multiplicative(Canonical::new()), factors, |_, e| {
        e as u64 + 1
    })
}

/// The sum of divisors `σ(n)` from a factorization of `n`.
///
/// # Definition
/// `σ(n) = Σ_{m: m|n} m`. It is multiplicative with `σ(p^e) = 1 + p + ... + p^e`.
///
/// # Contract
/// `factors` is a factorization, and `σ(n)` fits in `u64`.
///
/// # Complexity
/// - Time: O(log n)
/// - Space: O(1)
pub fn sum_divisors(factors: &[(u64, u32)]) -> u64 {
    multiplicative(&Multiplicative(Canonical::new()), factors, |p, e| {
        (0..=e).fold((0, 1), |(s, q), _| (s + q, q * p)).0
    })
}
