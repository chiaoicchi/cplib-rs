use crate::algebra::{Group, Monoid, Semigroup};
use crate::arithmetic::lpf::Lpf;

/// The symmetric group `S_n`.
///
/// # Definition
/// `S_n` is the group of the permutations of `[0, n)` under composition, where a permutation `p` is
/// stored as the vector with `p[i]` the image of `i`, and `op(p, q)` sends `i` to `q[p[i]]`, `id()`
/// is the identity, and `inv(p)` is the inverse.
///
/// # Contract
/// The values are permutations of `[0, n)`.
///
/// # Complexity
/// - Time: O(n) for `op`, `id` and `inv`.
/// - Space: O(n)
#[derive(Clone, Copy)]
pub struct Symmetric(pub usize);
impl Symmetric {
    /// The degree `n`.
    pub fn degree(&self) -> usize {
        self.0
    }
}
impl Semigroup for Symmetric {
    type Value = Vec<usize>;
    fn op(&self, p: &Self::Value, q: &Self::Value) -> Self::Value {
        p.iter().map(|&i| q[i]).collect()
    }
}
impl Monoid for Symmetric {
    fn id(&self) -> Self::Value {
        (0..self.0).collect()
    }
}
impl Group for Symmetric {
    fn inv(&self, p: &Self::Value) -> Self::Value {
        let mut q = vec![0; self.0];
        for (i, a) in p.iter().enumerate() {
            q[*a] = i;
        }
        q
    }
}

/// The cycles of `p`.
///
/// # Definition
/// Each cycle is `(i, p[i], p[p[i]], ...)` starting from its least element, and the cycles are in
/// increasing order of that element. Fixed points are cycles of length `1`.
///
/// # Contract
/// `p` is a permutation of `[0, n)`.
///
/// # Complexity
/// - Time: O(n)
/// - Space: O(n)
pub fn sym_cycles(p: &[usize]) -> Vec<Vec<usize>> {
    let n = p.len();
    let mut seen = vec![false; n];
    let mut cycles = Vec::new();
    for mut i in 0..n {
        if seen[i] {
            continue;
        }
        let mut cycle = Vec::new();
        while !seen[i] {
            seen[i] = true;
            cycle.push(i);
            i = p[i];
        }
        cycles.push(cycle);
    }
    cycles
}

/// Whether `p` is even.
///
/// # Definition
/// The sign of `p` is `(-1)^(n-c)`, where `c` is the number of cycles, and `p` is even iff the sign
/// is `1`.
///
///# Contract
/// `p` is a permutation of `[0, n)`.
//
/// # Complexity
/// - Time: O(n)
/// - Space: O(n)
pub fn sym_sign(p: &[usize]) -> bool {
    (p.len() - sym_cycles(p).len()) & 1 == 0
}

/// The order of `p`, as its factorization.
///
/// # Definition
/// The least `k >= 1` with `p^k = id`, which is the least common multiple of the cycle lengths, as
/// the pairs `(q, e)` of its prime factorization in increasing order of `q`.
///
/// # Contract
/// `p` is a permutation of `[0, n)`.
///
/// # Complexity
/// - Time: O(n log log n)
/// - Space: O(n)
pub fn sym_order(p: &[usize]) -> Vec<(usize, u32)> {
    let n = p.len();
    let lpf = Lpf::new(n);
    let mut exponent = std::collections::BTreeMap::new();
    let mut seen = vec![false; n];
    for mut i in 0..n {
        if seen[i] {
            continue;
        }
        let mut cnt = 0;
        while !seen[i] {
            seen[i] = true;
            cnt += 1;
            i = p[i];
        }
        for (q, e) in lpf.prime_factors(cnt) {
            let f = exponent.entry(q).or_insert(0);
            *f = (*f).max(e);
        }
    }
    exponent.into_iter().collect()
}

/// The power `p^k`.
///
/// # Definition
/// `p^0 = id` and `p^k = p^{k-1} p`.
///
/// # Contract
/// `p` is a permutation of `[0, n)`.
///
/// # Complexity
/// - Time: O(n)
/// - Space: O(n)
pub fn sym_pow(p: &[usize], k: u64) -> Vec<usize> {
    let n = p.len();
    let mut q = vec![0; n];
    for cycle in sym_cycles(p) {
        let l = cycle.len();
        let shift = (k % l as u64) as usize;
        for (i, &x) in cycle.iter().enumerate() {
            q[x] = cycle[(i + shift) % l];
        }
    }
    q
}
