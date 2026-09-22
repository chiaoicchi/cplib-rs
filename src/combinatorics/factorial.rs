use crate::algebra::Field;

/// The factorials `0!, 1!, ..., n!` and their inverses in `R`.
///
/// # Definition
/// For `i` in `[0, n]`, `i!` denotes its image in `R`.
///
/// # Invariants
/// - `fact[i] = i!` and `inv_fact[i] = (i!)^{-1}` for `i` in `[0, n]`.
///
/// # Complexity
/// - Space: O(n)
pub struct Factorial<R: Field> {
    ring: R,
    fact: Box<[R::Value]>,
    inv_fact: Box<[R::Value]>,
}

impl<R: Field<Value: Clone>> Factorial<R> {
    /// The factorials `0!, 1!, ..., n!` and their inverses in `R`.
    ///
    /// # Contract
    /// `1, ..., n` are invertible in `R`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    pub fn new(ring: R, n: usize) -> Self {
        let one = ring.one();
        let minus_one = ring.neg(&one);
        let mut fact = Vec::with_capacity(n + 1);
        fact.push(ring.one());
        let mut k = ring.zero();
        for i in 1..=n {
            k = ring.add(&k, &one);
            fact.push(ring.mul(&fact[i - 1], &k));
        }
        let mut inv_fact: Vec<R::Value> = Vec::with_capacity(n + 1);
        let mut inv = ring.inv(&fact[n]);
        for _ in 0..n {
            let next = ring.mul(&inv, &k);
            inv_fact.push(inv);
            inv = next;
            k = ring.add(&k, &minus_one);
        }
        inv_fact.push(inv);
        inv_fact.reverse();
        Self {
            ring,
            fact: fact.into(),
            inv_fact: inv_fact.into(),
        }
    }

    /// The factorial `i!`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i > n`.
    pub fn factorial(&self, i: usize) -> R::Value {
        assert!(
            i < self.fact.len(),
            "index out of bounds: i={i}, n={}",
            self.fact.len() - 1,
        );
        self.fact[i].clone()
    }

    /// The inverse `(i!)^{-1}` of the factorial.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i > n`.
    pub fn inv_factorial(&self, i: usize) -> R::Value {
        assert!(
            i < self.fact.len(),
            "index out of bounds: i={i}, n={}",
            self.fact.len() - 1,
        );
        self.inv_fact[i].clone()
    }

    /// The inverse `i^{-1}`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i = 0` or `i > n`.
    pub fn inv(&self, i: usize) -> R::Value {
        assert!(i != 0, "zero has no inverse");
        assert!(
            i < self.fact.len(),
            "index out of bounds: i={i}, n={}",
            self.fact.len() - 1,
        );
        self.ring.mul(&self.fact[i - 1], &self.inv_fact[i])
    }

    /// The binomial coefficient of `a` and `b`.
    ///
    /// # Definition
    /// `[x^b](1 + x)^a`, which is `a! / (b! (a - b)!)` for `b <= a` and `0` otherwise.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `a > n`.
    pub fn binomial(&self, a: usize, b: usize) -> R::Value {
        assert!(
            a < self.fact.len(),
            "index out of bounds: a={a}, n={}",
            self.fact.len() - 1,
        );
        if b > a {
            self.ring.zero()
        } else {
            self.ring.mul(
                &self.ring.mul(&self.fact[a], &self.inv_fact[a - b]),
                &self.inv_fact[b],
            )
        }
    }

    /// The falling factorial of `a` of length `b`.
    ///
    /// # Definition
    /// `a (a - 1) ... (a - b + 1)`, the product of `b` factors. It is `0` for `b > a`, as one of
    /// the factors is `0`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `a > n`.
    pub fn falling_factorial(&self, a: usize, b: usize) -> R::Value {
        assert!(
            a < self.fact.len(),
            "index out of bounds: a={a}, n={}",
            self.fact.len() - 1,
        );
        if b > a {
            self.ring.zero()
        } else {
            self.ring.mul(&self.fact[a], &self.inv_fact[a - b])
        }
    }

    /// The multinomial coefficient of `ks`.
    ///
    /// # Definition
    /// For `ks = (k_1, ..., k_m)` and `s = k_1 + ... + k_m`, the coefficient
    /// `[x_1^{k_1} ... x_m^{k_m}] (x_1 + ... + x_m)^s`, which is `s! / (k_1! ... k_m!)`.
    ///
    /// # Complexity
    /// - Time: O(m)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `s > n`.
    pub fn multinomial(&self, ks: &[usize]) -> R::Value {
        let s: usize = ks.iter().sum();
        assert!(
            s < self.fact.len(),
            "index out of bounds: ks={:?}, n={}",
            ks,
            self.fact.len() - 1,
        );
        ks.iter().fold(self.fact[s].clone(), |acc, a| {
            self.ring.mul(&acc, &self.inv_fact[*a])
        })
    }
}
