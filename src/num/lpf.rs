/// The least prime factor of every integer in `[0, n]`, by the linear sieve.
///
/// # Definition
/// For `m >= 2`, `lpf(m)` is the least prime dividing `m`; it exists by unique factorization, and
/// `lpf(m) == m` exactly when `m` is prime. Dividing `lpf(m)` out repeatedly strictly decreases
/// `m`, so the table determines the factorization of every `m` in `[1, n]` in `Ω(m)` divisions,
/// where `Ω(m)` is the number of prime factors of `m` with multiplicity and `ω(m)` the number of
/// distinct ones.
///
/// # Invariants
/// - `value[m]` is the least prime dividing `m` for `2 <= m <= n`; `value[0] = 0` and
///   `value[1] = 1`, neither of which is a prime factor.
///
/// # Complexity
/// - Space: O(n)
pub struct Lpf {
    value: Box<[u32]>,
}

impl Lpf {
    /// Constructs the table for `[0, n]`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `n` is greater than or equal to `2^32`.
    pub fn new(n: usize) -> Self {
        assert!(n < 1 << 32, "n must be less than 2^32: n={n}");
        let mut value = vec![0; n + 1];
        if n > 1 {
            value[1] = 1;
        }
        let mut ps = Vec::new();
        for i in 2..=n {
            if value[i] == 0 {
                value[i] = i as u32;
                ps.push(i);
            }
            for &p in &ps {
                if p as u32 > value[i] || i * p > n {
                    break;
                }
                value[i * p] = p as u32;
            }
        }
        Self {
            value: value.into(),
        }
    }

    /// The largest `m` the table answers for.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn bound(&self) -> usize {
        self.value.len() - 1
    }

    /// The least prime dividing `m`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `m < 2` or `m > bound()`.
    pub fn lpf(&self, m: usize) -> usize {
        assert!(m >= 2, "m must be greater than 1: m={m}");
        assert!(
            m <= self.bound(),
            "m must be less than or equal to bound: m={m}, bound={}",
            self.bound()
        );
        self.value[m] as usize
    }

    /// Returns `true` if `m` is prime.
    ///
    /// # Definition
    /// `m` is prime iff `lpf(m) = m`; `0` and `1` are not prime.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `m > bound()`
    pub fn is_prime(&self, m: usize) -> bool {
        assert!(
            m <= self.bound(),
            "m must be less than or equal to bound: m={m}, bound={}",
            self.bound()
        );
        m >= 2 && self.value[m] as usize == m
    }

    /// The prime factorization of `m`, in increasing order of the prime.
    ///
    /// # Definition
    /// Yields the pairs `(p, e)` with `p` prime, `e >= 1` and `m = Π p^e`, in increasing order of
    /// `p`; `m = 1` yields nothing. The primes come out in order because dividing out the least
    /// prime factor leaves a cofactor all of whose prime factors are at least as large. No pair is
    /// stored, so an iteration allocates nothing.
    ///
    /// # Complexity
    /// - Time: O(Ω(m)) in total, O(1) amortized per item
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `m == 0` or `m > bound()`.
    pub fn prime_factors(&self, mut m: usize) -> impl Iterator<Item = (usize, u32)> + '_ {
        assert!(m > 0, "m must be greater than 0: m={m}");
        assert!(
            m <= self.bound(),
            "m must be less than or equal to bound: m={m}, bound={}",
            self.bound()
        );
        std::iter::from_fn(move || {
            if m == 1 {
                return None;
            }
            let p = self.value[m] as usize;
            let mut e = 0;
            while m % p == 0 {
                m /= p;
                e += 1;
            }
            Some((p, e))
        })
    }

    /// The prime factorization of `m`, collected.
    ///
    /// # Complexity
    /// - Time: O(Ω(m))
    /// - Space: O(ω(m))
    ///
    /// # Panics
    /// Panics if `m == 0` or `m > bound()`.
    pub fn factorize(&self, m: usize) -> Vec<(usize, u32)> {
        self.prime_factors(m).collect()
    }

    /// The primes in `[2, bound()]` in increasing order.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(1)
    pub fn primes(&self) -> impl Iterator<Item = usize> + '_ {
        (2..=self.bound()).filter(|i| self.value[*i] == *i as u32)
    }

    /// The divisors of `m`, in no particular order.
    ///
    /// # Definition
    /// A divisor of `m = Π p^e` is `Π p^f` with `0 <= f <= e`, so the divisors are obtained by
    /// multiplying the list by each `p^f` in turn; there are `d(m) = Π (e + 1)` of them. They come
    /// out in the order that construction produces, which is not increasing.
    ///
    /// # Complexity
    /// - Time: O(Ω(m) + d(m))
    /// - Space: O(d(m))
    ///
    /// # Panics
    /// Panics if `m = 0` or `m > bound()`.
    pub fn divisors(&self, m: usize) -> Vec<usize> {
        let mut ds = vec![1];
        for (p, e) in self.prime_factors(m) {
            let len = ds.len();
            let mut q = 1;
            for _ in 0..e {
                q *= p;
                for i in 0..len {
                    ds.push(ds[i] * q);
                }
            }
        }
        ds
    }

    /// Returns `true` if `m` is squarefree.
    ///
    /// # Definition
    /// `m` is squarefree if no `p^2` divides it, that is if every exponent of its factorization is
    /// `1`; equivalently `m = radical(m)`, equivalently `μ(m) != 0`. `1` is squarefree.
    ///
    /// # Complexity
    /// - Time: O(Ω(m))
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `m = 0` or `m > bound()`.
    pub fn is_squarefree(&self, m: usize) -> bool {
        self.prime_factors(m).all(|x| x.1 == 1)
    }

    /// The radical of `m`.
    ///
    /// # Definition
    /// `rad(m) = Π_{p | m} p`, the product of the distinct primes dividing `m`; `rad(1) = 1`. It is
    /// the largest squarefree divisor of `m`, and `rad(a) = rad(b)` iff `a` and `b` have the same
    /// prime divisors.
    ///
    /// # Complexity
    /// - Time: O(Ω(m))
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `m = 0` or `m > bound()`.
    pub fn radical(&self, m: usize) -> usize {
        self.prime_factors(m).map(|(p, _)| p).product()
    }

    /// The number of distinct primes dividing `m`.
    ///
    /// # Definition
    /// `ω(m)`, the number of pairs of the factorization; `ω(1) = 0`.
    ///
    /// # Complexity
    /// - Time: O(Ω(m))
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `m = 0` or `m > bound()`.
    pub fn omega(&self, m: usize) -> u32 {
        self.prime_factors(m).count() as u32
    }

    /// The number of primes dividing `m` with multiplicity.
    ///
    /// # Definition
    /// `Ω(m) = Σ e` over the factorization; `Ω(1) = 0`. It is the length of a maximal chain of
    /// proper divisors of `m`, and `(-1)^Ω` is the Liouville function.
    ///
    /// # Complexity
    /// - Time: O(Ω(m))
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `m = 0` or `m > bound()`.
    pub fn big_omega(&self, m: usize) -> u32 {
        self.prime_factors(m).map(|(_, e)| e).sum()
    }
}
