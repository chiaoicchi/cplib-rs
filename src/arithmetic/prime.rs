/// Whether `n` is prime.
///
/// # Definition
/// `n` is prime if `n >= 2` and its only positive divisors are `1` and `n`.
///
/// # Complexity
/// - Time: O(√n)
/// - Space: O(1)
pub const fn is_prime(n: u32) -> bool {
    if n < 2 {
        false
    } else if n & 1 == 0 {
        n == 2
    } else {
        let mut d = 3;
        while d <= n / d {
            if n % d == 0 {
                return false;
            }
            d += 2;
        }
        true
    }
}

/// The primes in `[2, n]`, in increasing order.
///
/// # Definition
/// The `p` in `[2, n]` whose only positive divisors are `1` and `p`.
///
/// # Complexity
/// - Time: O(n log log n)
/// - Space: O(n)
pub fn primes(n: usize) -> Vec<usize> {
    if n < 2 {
        return Vec::new();
    }
    let half = (n - 1) >> 1;
    let mut composite = vec![true; half + 1];
    let mut ps = vec![2];
    for i in 1..=half {
        if !composite[i] {
            continue;
        }
        let p = (i << 1) + 1;
        ps.push(p);
        let mut j = (p * p - 1) >> 1;
        while j <= half {
            composite[j] = true;
            j += p;
        }
    }
    ps
}

/// The prime factorization of `n`.
///
/// # Definition
/// The pairs `(p, e)` with `p` prime, `e >= 1` and `n = Π p^e`, in increasing order of `p`. For
/// `n = 1` it is empty.
///
/// # Complexity
/// - Time: O(√n)
/// - Space: O(log n)
///
/// # Panics
/// Panics if `n == 0`.
pub fn factorize(mut n: u64) -> Vec<(u64, u32)> {
    assert!(n >= 1, "n must be positive");
    let mut factors = Vec::new();
    let mut p = 2;
    while p <= n / p {
        if n % p == 0 {
            let mut e = 0;
            while n % p == 0 {
                n /= p;
                e += 1;
            }
            factors.push((p, e));
        }
        p += if p == 2 { 1 } else { 2 };
    }
    if n > 1 {
        factors.push((n, 1));
    }
    factors
}
