use crate::algebra::{Inv, One, Zero};

/// The table of factorials `0!, 1!, ..., n!` and their inverses, represented in `T`.
///
/// # Definition
/// `factorial(i)` is the image in `T` of `i!`, and `inv_factorial(i)` is its inverse, for `i` in
/// `[0, n]`. The binomial coefficient is `binomial(a, b) = [x^b](1 + x)^a`, which is
/// `a! / (b! (a - b)!)` for `b <= a` and `0` otherwise.
///
/// # Invariants
/// `fact[i]` and `inv_fact[i]` are the images of `i!` and `(i!)^{-1}`, and `fact[i] * inv_fact[i]`
/// is `one` for `i` in `[0, n]`.
///
/// # Complexity
/// - Space: O(n)
pub struct Factorial<T> {
    fact: Box<[T]>,
    inv_fact: Box<[T]>,
}

impl<
    T: Clone
        + Zero
        + One
        + Inv<Output = T>
        + std::ops::Add<Output = T>
        + std::ops::Sub<Output = T>
        + std::ops::Mul<Output = T>,
> Factorial<T>
{
    /// Constructs the table for `i` in `[0, n]`.
    ///
    /// `n!` must be invertible in `T`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    pub fn new(n: usize) -> Self {
        let mut fact = Vec::with_capacity(n + 1);
        let mut inv_fact = Vec::with_capacity(n + 1);
        fact.push(T::one());
        let mut t = T::one();
        let mut x = T::one();
        for _ in 1..=n {
            x = x * t.clone();
            fact.push(x.clone());
            t = t + T::one();
        }
        x = x.inv();
        inv_fact.push(x.clone());
        for _ in (1..=n).rev() {
            t = t - T::one();
            x = x * t.clone();
            inv_fact.push(x.clone());
        }
        inv_fact.reverse();
        Self {
            fact: fact.into(),
            inv_fact: inv_fact.into(),
        }
    }

    /// Returns `i!`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i > n`.
    pub fn factorial(&self, i: usize) -> T {
        assert!(
            i < self.fact.len(),
            "index out of bounds: i={i}, n={}",
            self.fact.len() - 1,
        );
        self.fact[i].clone()
    }

    /// Returns `(i!)^{-1}`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i > n`.
    pub fn inv_factorial(&self, i: usize) -> T {
        assert!(
            i < self.fact.len(),
            "index out of bounds: i={i}, n={}",
            self.fact.len() - 1,
        );
        self.inv_fact[i].clone()
    }

    /// Returns `i^{-1}`, as `(i - 1)! (i!)^{-1}`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i == 0` or `i > n`.
    pub fn inv(&self, i: usize) -> T {
        assert!(i != 0, "zero has no inverse");
        assert!(
            i < self.fact.len(),
            "index out of bounds: i={i}, n={}",
            self.fact.len() - 1,
        );
        self.fact[i - 1].clone() * self.inv_fact[i].clone()
    }

    /// Returns `binomial(a, b) = [x^b](1 + x)^a`, which is
    /// `a (a - 1) ... (a - b + 1) / 1 2 ... b = a! / (b! (a - b)!)` for
    /// `b <= a` and `0` otherwise.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `a > n`.
    pub fn binomial(&self, a: usize, b: usize) -> T {
        assert!(
            a < self.fact.len(),
            "index out of bounds: a={a}, n={}",
            self.fact.len() - 1,
        );
        if b > a {
            T::zero()
        } else {
            self.fact[a].clone() * self.inv_fact[a - b].clone() * self.inv_fact[b].clone()
        }
    }

    /// Returns the falling factorial `b! [x^b] (1 + x)^a`, which is
    /// `a (a - 1) ... (a - b + 1) = a! / (a - b)!` for `b <= a` and `0` otherwise.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `a > n`.
    pub fn falling_factorial(&self, a: usize, b: usize) -> T {
        assert!(
            a < self.fact.len(),
            "index out of bounds: a={a}, n={}",
            self.fact.len() - 1,
        );
        if b > a {
            T::zero()
        } else {
            self.fact[a].clone() * self.inv_fact[a - b].clone()
        }
    }

    /// Returns the multinomial coefficient
    /// `[x_1^{k_1} ... x_m^{k_m}] (x_1 + ... + x_m)^s` with `s = k_1 + ... + k_m`, which is
    /// `s! / (k_1! ... k_m!)`.
    ///
    /// # Complexity
    /// - Time: O(m)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `k_1 + ... + k_m > n`.
    pub fn multinomial(&self, ks: &[usize]) -> T {
        let s: usize = ks.iter().sum();
        assert!(
            s < self.fact.len(),
            "index out of bounds: ks={:?}, n={}",
            ks,
            self.fact.len() - 1,
        );
        ks.iter().fold(self.fact[s].clone(), |acc, a| {
            acc * self.inv_fact[*a].clone()
        })
    }
}
