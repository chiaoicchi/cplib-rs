/// Returns `true` if `n` is prime.
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

/// The primes below `n` by the sieve of Eratosthenes.
///
/// # Complexity
/// - Time: O(n log log n)
/// - Space: O(n)
pub fn primes(n: usize) -> Vec<usize> {
    let mut is_prime = vec![true; n + 1];
    let mut ps = Vec::new();
    for p in 2..=n {
        if !is_prime[p] {
            continue;
        }
        ps.push(p);
        let mut m = p * p;
        while m <= n {
            is_prime[m] = false;
            m += p;
        }
    }
    ps
}

/// The prime factorization of `n` by trial division.
///
/// # Definition
/// Returns the pairs `(p, e)` with `p` prime, `e >= 1` and `n = Π p^e`, in increasing order
/// of `p`. By unique factorization this representation is unique; `n = 1` gives the empty product.
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

/// The smallest prime factor of every integer below `n` by the linear sieve.
///
/// # Definition
/// `spf[m]` is the least prime dividing `m` for `m >= 2`; `spf[m] == m` iff `m` is prime.
/// `spf[0] = 0` and `spf[1] = 1`; neither is a prime factor.
///
/// # Complexity
/// - Time: O(n)
/// - Space: O(n)
pub fn spf(n: usize) -> Vec<usize> {
    let mut spf = vec![0; n + 1];
    if n > 1 {
        spf[1] = 1;
    }
    let mut ps = Vec::new();
    for i in 2..=n {
        if spf[i] == 0 {
            spf[i] = i;
            ps.push(i);
        }
        for &p in &ps {
            if p > spf[i] || i * p > n {
                break;
            }
            spf[i * p] = p;
        }
    }
    spf
}
