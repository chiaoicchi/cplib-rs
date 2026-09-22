/// The least prime factors of the integers in `[0, n]`.
///
/// # Definition
/// For `m >= 2`, `lpf(m)` is the least prime dividing `m`. In the complexities below, `Ω(m)` and
/// `ω(m)` are the numbers of prime factors of `m` with and without multiplicity, and `d(m)` is the
/// number of divisors of `m`.
///
/// # Invariants
/// - `value[m] = lpf(m)` for `2 <= m <= n`, `value[0] = 0` and `value[1] = 1` if `n >= 1`.
///
/// # Complexity
/// - Space: O(n)
pub struct Lpf {
    value: Box<[u32]>,
}

impl Lpf {
    /// The least prime factors of the integers in `[0, n]`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    ///
    /// # Panics
    /// Panics if `n >= 2^32`.
    pub fn new(n: usize) -> Self {
        assert!(n < 1 << 32, "n must be less than 2^32: n={n}");
        let mut value = vec![0; n + 1];
        if n > 0 {
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

    /// The bound `n` of the table.
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

    /// Whether `m` is prime.
    ///
    /// # Definition
    /// `m` is prime if `m >= 2` and its only positive divisors are `1` and `m`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `m > n`.
    pub fn is_prime(&self, m: usize) -> bool {
        assert!(
            m <= self.bound(),
            "m must be less than or equal to bound: m={m}, bound={}",
            self.bound()
        );
        m >= 2 && self.value[m] as usize == m
    }

    /// The prime factorization of `m`, as an iterator.
    ///
    /// # Definition
    /// The pairs `(p, e)` with `p` prime, `e >= 1` and `m = Π p^e`, in increasing order of `p`. For
    /// `m = 1` it is empty.
    ///
    /// # Complexity
    /// - Time: O(Ω(m)) in total
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `m = 0` or `m > n`.
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

    /// The prime factorization of `m`, as a vector.
    ///
    /// # Definition
    /// The pairs of [`Lpf::prime_factors`], in the same order.
    ///
    /// # Complexity
    /// - Time: O(Ω(m))
    /// - Space: O(ω(m))
    ///
    /// # Panics
    /// Panics if `m = 0` or `m > n`.
    pub fn factorize(&self, m: usize) -> Vec<(usize, u32)> {
        self.prime_factors(m).collect()
    }

    /// The primes in `[2, n]` in increasing order.
    ///
    /// # Definition
    /// The `p` in `[2, n]` whose only positive divisors are `1` and `p`.
    ///
    /// # Complexity
    /// - Time: O(n) in total
    /// - Space: O(1)
    pub fn primes(&self) -> impl Iterator<Item = usize> + '_ {
        (2..=self.bound()).filter(|i| self.value[*i] == *i as u32)
    }

    /// The divisors of `m`.
    ///
    /// # Definition
    /// Every positive divisor of `m` exactly once, in an unspecified order.
    ///
    /// # Complexity
    /// - Time: O(Ω(m) + d(m))
    /// - Space: O(d(m))
    ///
    /// # Panics
    /// Panics if `m = 0` or `m > n`.
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

    /// Whether `m` is squarefree.
    ///
    /// # Definition
    /// `m` is squarefree if `p^2` divides `m` for no prime `p`. `1` is squarefree.
    ///
    /// # Complexity
    /// - Time: O(Ω(m))
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `m = 0` or `m > n`.
    pub fn is_squarefree(&self, m: usize) -> bool {
        self.prime_factors(m).all(|x| x.1 == 1)
    }

    /// The radical of `m`.
    ///
    /// # Definition
    /// `rad(m) = Π_{p | m} p`, the product of the distinct primes dividing `m`. `rad(1) = 1`.
    ///
    /// # Complexity
    /// - Time: O(Ω(m))
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `m = 0` or `m > n`.
    pub fn radical(&self, m: usize) -> usize {
        self.prime_factors(m).map(|(p, _)| p).product()
    }

    /// The number `ω(m)` of distinct primes dividing `m`.
    ///
    /// # Definition
    /// `ω(m)`, the number of pairs of the factorization; `ω(1) = 0`.
    ///
    /// # Complexity
    /// - Time: O(Ω(m))
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `m = 0` or `m > n`.
    pub fn omega(&self, m: usize) -> u32 {
        self.prime_factors(m).count() as u32
    }

    /// The number `Ω(m)` of primes dividing `m`, counted with multiplicity.
    ///
    /// # Definition
    /// `Ω(m) = Σ e` for `m = Π p^e`. `Ω(1) = 0`.
    ///
    /// # Complexity
    /// - Time: O(Ω(m))
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `m = 0` or `m > n`.
    pub fn big_omega(&self, m: usize) -> u32 {
        self.prime_factors(m).map(|(_, e)| e).sum()
    }
}
