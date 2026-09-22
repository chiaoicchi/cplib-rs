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
    let mut composite = vec![false; half + 1];
    let mut ps = vec![2];
    for i in 1..=half {
        if composite[i] {
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

/// Whether `2^k - 1` is prime.
///
/// # Complexity
/// - Time: O(k)
/// - Space: O(1)
///
/// # Panics
/// Panics if `k > 64`.
pub const fn is_mersenne_exponent(k: u32) -> bool {
    assert!(k <= 64, "k must be at most 64");
    if k == 2 {
        return true;
    }
    if k < 2 || !is_prime(k) {
        return false;
    }
    let m = (1u128 << k) - 1;
    let mut s = 4;
    let mut i = 0;
    while i < k - 2 {
        s = (s * s + m - 2) % m;
        i += 1;
    }
    s == 0
}

/// The power `base^exp` modulo `p`, in `[0, p)`.
///
/// # Contract
/// `2 <= p <= 2^32`.
///
/// # Complexity
/// - Time: O(log exp)
/// - Space: O(1)
const fn pow_mod(mut base: u64, mut exp: u64, p: u64) -> u64 {
    let mut x = 1;
    base %= p;
    while exp > 0 {
        if exp & 1 == 1 {
            x = x * base % p;
        }
        base = base * base % p;
        exp >>= 1;
    }
    x
}

/// The least primitive root modulo `p`.
///
/// # Definition
/// The least generator `g` of the multiplicative group `(Z/pZ)^*`, as an integer in `[1, p)`. For
/// `p = 2` it is `1`.
///
/// Equivalently, `g` is the least integer in `[1, p)` whose order modulo `p` is `p - 1`.
///
/// # Complexity
/// - Time: O(√p + g k log p), where `k` is the number of prime factors of
///   `p - 1`
/// - Space: O(1)
///
/// # Panics
/// Panics if `p` is not prime.
pub const fn primitive_root(p: u32) -> u32 {
    assert!(is_prime(p), "p must be prime");
    if p == 2 {
        return 1;
    }
    let mut factors = [0u32; 32];
    let mut k = 0;
    let mut m = p - 1;
    let mut d = 2;
    while d as u64 * d as u64 <= m as u64 {
        if m % d == 0 {
            factors[k] = d;
            k += 1;
            while m % d == 0 {
                m /= d;
            }
        }
        d += 1;
    }
    if m > 1 {
        factors[k] = m;
        k += 1;
    }
    let mut g = 2;
    loop {
        let mut f = true;
        let mut i = 0;
        while i < k {
            if pow_mod(g as u64, ((p - 1) / factors[i]) as u64, p as u64) == 1 {
                f = false;
                break;
            }
            i += 1;
        }
        if f {
            return g;
        }
        g += 1;
    }
}
