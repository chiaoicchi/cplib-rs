use crate::range::to_half_open;

/// A subset `S` of `[0, u)`, as a 64-ary tree of bit words.
///
/// # Definition
/// `u` us the bound of `S`, `n = |S|`.
///
/// # Invariants
/// - `d` is the least `d >= 1` with `64^d >= u`, and word[level[h]..level[h + 1]]` is the level
///   `h` for `h` in `[0, d)`, of `ceil(u / 64^{h+1})` words.
/// - Bit `j` of `word[level[0] + i]` is set iff `64i + j` is in `S`, and bit `j` if
///   `word[level[h + 1] + i]` is set iff `word[level[h] + 64i + j]` is nonzero.
/// - `len = n`.
///
/// # Complexity
/// - Space: O(u / 64)
pub struct IntSet {
    word: Box<[u64]>,
    level: Box<[usize]>,
    len: usize,
    bound: usize,
}

impl IntSet {
    /// The empty subset of `[0, u)`
    ///
    /// # Complexity
    /// - Time: O(u / 64)
    /// - Space: O(u / 64)
    pub fn new(u: usize) -> Self {
        let mut level = vec![0];
        let mut p = u;
        loop {
            let w = p.div_ceil(64);
            level.push(level[level.len() - 1] + w);
            p = w;
            if p <= 1 {
                break;
            }
        }
        Self {
            word: vec![0; level[level.len() - 1]].into(),
            level: level.into(),
            len: 0,
            bound: u,
        }
    }

    /// The subset `[0, u)` itself.
    ///
    /// # Complexity
    /// - Time: O(u / 64)
    /// - Space: O(u / 64)
    pub fn full(u: usize) -> Self {
        let mut x = Self::new(u);
        let mut p = u;
        for h in 0..x.level.len() - 1 {
            let (lo, hi) = (x.level[h], x.level[h + 1]);
            for k in lo..hi {
                x.word[k] = !0;
            }
            if p & 63 != 0 {
                x.word[lo + (p >> 6)] = (1 << (p & 63)) - 1;
            }
            p = p.div_ceil(64);
        }
        x.len = u;
        x
    }

    /// The subset of `[0, u)` of the elements of `xs`.
    ///
    /// # Complexity
    /// - Time: O(u / 64 + xs.len())
    /// - Space: O(u / 64)
    ///
    /// # Panics
    /// Panics if `xs[i] >= u` for some `i`.
    pub fn from_slice(u: usize, xs: &[usize]) -> Self {
        let mut x = Self::new(u);
        for &v in xs {
            assert!(v < u, "out of bounds: x={v}, u={u}");
            x.word[v >> 6] |= 1 << (v & 63);
        }
        x.len = (0..x.level[1])
            .map(|k| x.word[k].count_ones() as usize)
            .sum();
        for h in 1..x.level.len() - 1 {
            for k in x.level[h - 1]..x.level[h] {
                if x.word[k] != 0 {
                    let p = k - x.level[h - 1];
                    x.word[x.level[h] + (p >> 6)] |= 1 << (p & 63);
                }
            }
        }
        x
    }

    /// Inserts `x` into `S`, and returns whether `S` changed.
    ///
    /// # Complexity
    /// - Time: O(log_64 u)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x >= u`.
    pub fn insert(&mut self, x: usize) -> bool {
        assert!(x < self.bound, "out of bounds: x={x}, u={}", self.bound);
        if self.word[x >> 6] >> (x & 63) & 1 == 1 {
            return false;
        }
        self.len += 1;
        let mut p = x;
        for h in 0..self.level.len() - 1 {
            let k = self.level[h] + (p >> 6);
            let empty = self.word[k] == 0;
            self.word[k] |= 1 << (p & 63);
            if !empty {
                break;
            }
            p >>= 6;
        }
        true
    }

    /// Removes `x` from `S`, and returns whether `S` changed.
    ///
    /// # Complexity
    /// - Time: O(log_64 u)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x >= u`.
    pub fn remove(&mut self, x: usize) -> bool {
        assert!(x < self.bound, "out of bounds: x={x}, u={}", self.bound);
        if self.word[x >> 6] >> (x & 63) & 1 == 0 {
            return false;
        }
        self.len -= 1;
        let mut p = x;
        for h in 0..self.level.len() - 1 {
            let k = self.level[h] + (p >> 6);
            self.word[k] &= !(1 << (p & 63));
            if self.word[k] != 0 {
                break;
            }
            p >>= 6;
        }
        true
    }

    /// Whether `x` is in `S`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x >= u`.
    pub fn contains(&self, x: usize) -> bool {
        assert!(x < self.bound, "out of bounds: x={x}, u={}", self.bound);
        self.word[x >> 6] >> (x & 63) & 1 == 1
    }

    /// The least element of `S` at least `x`, or `None` if there is none.
    ///
    /// # Complexity
    /// - Time: O(log_64 u)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x > u`.
    pub fn min_ge(&self, x: usize) -> Option<usize> {
        assert!(x <= self.bound, "out of bounds: x={x}, u={}", self.bound);
        if x == self.bound {
            return None;
        }
        let mut p = x;
        for h in 0..self.level.len() - 1 {
            if self.level[h] + (p >> 6) >= self.level[h + 1] {
                break;
            }
            let w = self.word[self.level[h] + (p >> 6)] >> (p & 63);
            if w != 0 {
                p += w.trailing_zeros() as usize;
                for g in (0..h).rev() {
                    p = (p << 6) + self.word[self.level[g] + p].trailing_zeros() as usize;
                }
                return Some(p);
            }
            p = (p >> 6) + 1;
        }
        None
    }

    /// The greatest element of `S` at most `x`, or `None` if there is none.
    ///
    /// # Complexity
    /// - Time: O(log_64 u)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x >= u`.
    pub fn max_le(&self, x: usize) -> Option<usize> {
        assert!(x < self.bound, "out of bounds: x={x}, u={}", self.bound);
        let mut p = x;
        for h in 0..self.level.len() - 1 {
            let w = self.word[self.level[h] + (p >> 6)] << (63 - (p & 63));
            if w != 0 {
                p -= w.leading_zeros() as usize;
                for g in (0..h).rev() {
                    p = (p << 6) + 63 - self.word[self.level[g] + p].leading_zeros() as usize;
                }
                return Some(p);
            }
            if p >> 6 == 0 {
                break;
            }
            p = (p >> 6) - 1;
        }
        None
    }

    /// The least element of `S`, or `None` if `S` is empty.
    ///
    /// # Complexity
    /// - Time: O(log_64 u)
    /// - Space: O(1)
    pub fn min(&self) -> Option<usize> {
        (self.bound > 0).then(|| self.min_ge(0)).flatten()
    }

    /// The greatest element of `S`, or `None` if `S` is empty.
    ///
    /// # Complexity
    /// - Time: O(log_64 u)
    /// - Space: O(1)
    pub fn max(&self) -> Option<usize> {
        (self.bound > 0)
            .then(|| self.max_le(self.bound - 1))
            .flatten()
    }

    /// The elements of `S` in `range = [l, r)`, in increasing order.
    ///
    /// # Complexity
    /// - Time: O(log_64 u + k) to exhaust, where `k` is the number of elements yielded
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `l > r` or `r > u`.
    pub fn range(&self, range: impl std::ops::RangeBounds<usize>) -> impl Iterator<Item = usize> {
        let (l, r) = to_half_open(self.bound, range);
        let mut cur = 0u64;
        let mut base = 0;
        let mut start = l;
        std::iter::from_fn(move || {
            loop {
                if cur != 0 {
                    let x = base + cur.trailing_zeros() as usize;
                    cur &= cur - 1;
                    return (x < r).then_some(x);
                }
                if start >= r {
                    return None;
                }
                let x = self.min_ge(start)?;
                if x >= r {
                    return None;
                }
                base &= !63;
                cur = self.word[x >> 6] & !((1 << (x & 63)) - 1);
                start = base + 64;
            }
        })
    }

    /// The elements of `S` in increasing order.
    ///
    /// # Complexity
    /// - Time: O(k), where `k` is the number of elements in `S`
    /// - Space: O(1)
    pub fn iter(&self) -> impl Iterator<Item = usize> {
        self.range(..)
    }

    /// The size `n` of `S`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether `n = 0`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// The bound `u`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn bound(&self) -> usize {
        self.bound
    }
}
