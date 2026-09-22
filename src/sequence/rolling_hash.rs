use crate::algebra::{Monoid, Semigroup};
use crate::num::mersenne::Mersenne;
use crate::range::to_half_open;

/// A polynomial hash of a sequence with its length and base.
///
/// # Definition
/// For a base `b` in `Z/MZ` with `M = 2^K - 1`, the hash of a sequence `t` of length `m` is
/// `H(t) = t[0] b^(m-1) + t[1] b^(m-2) + ... + t[m-1]`, with each `t[i]` read in `Z/MZ`. Two hashes
/// are equal iff their values, lengths and bases are.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct SequenceHash<const K: u32> {
    value: Mersenne<K>,
    base: Mersenne<K>,
    len: usize,
}

impl<const K: u32> SequenceHash<K> {
    /// The value `H(t)`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn value(self) -> Mersenne<K> {
        self.value
    }

    /// The length `m` of `t`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(self) -> usize {
        self.len
    }

    /// Whether `m = 0`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(self) -> bool {
        self.len == 0
    }
}

/// The hash `H(t u)` of the concatenation, for `x = H(t)` and `y = H(u)`.
///
/// # Complexity
/// - Time: O(log y.len())
/// - Space: O(1)
///
/// # Panics
/// Panics if the bases of `x` and `y` differ.
pub fn sequence_concat<const K: u32>(x: SequenceHash<K>, y: SequenceHash<K>) -> SequenceHash<K> {
    assert!(x.base == y.base, "bases differ: x={}, y={}", x.base, y.base);
    SequenceHash {
        value: x.value * x.base.pow(y.len as u64) + y.value,
        base: x.base,
        len: x.len + y.len,
    }
}

/// The monoid of sequence hashes with the base `self.0` under [`sequence_concat`].
///
/// # Definition
/// `op(x, y) = concat(x, y)`, and `id()` is the hash of the empty sequence with the base `self.0`.
///
/// # Contract
/// The values have the base `self.0`.
#[derive(Clone, Copy)]
pub struct Concat<const K: u32>(pub Mersenne<K>);
impl<const K: u32> Semigroup for Concat<K> {
    type Value = SequenceHash<K>;
    fn op(&self, a: &SequenceHash<K>, b: &SequenceHash<K>) -> SequenceHash<K> {
        sequence_concat(*a, *b)
    }
}
impl<const K: u32> Monoid for Concat<K> {
    fn id(&self) -> SequenceHash<K> {
        SequenceHash {
            value: Mersenne::new(0),
            base: self.0,
            len: 0,
        }
    }
}

/// The hashes of the contiguous subsequences of a sequence.
///
/// # Definition
/// For a sequence `s` of length `n` and a base `b`, gives `H(s[l..r))` for `0 <= l <= r <= n`, with
/// each `s[i]` read in `Z/MZ` through `u64`.
///
/// # Invariants
/// - `prefix[i] = H(s[0..i))` and `power[i] = b^i` for `i` in `[0, n]`, and `base = b`.
///
/// # Complexity
/// - Space: O(n)
pub struct RollingHash<const K: u32> {
    prefix: Box<[Mersenne<K>]>,
    power: Box<[Mersenne<K>]>,
    base: Mersenne<K>,
}

impl<const K: u32> RollingHash<K> {
    /// The rolling hash of `s` with the base `base`.
    ///
    /// # Complexity
    /// - Time: O(n)
    /// - Space: O(n)
    pub fn from_slice<T: Copy + Into<u64>>(base: Mersenne<K>, s: &[T]) -> Self {
        let n = s.len();
        let mut prefix = Vec::with_capacity(n + 1);
        let mut power = Vec::with_capacity(n + 1);
        prefix.push(Mersenne::new(0));
        power.push(Mersenne::new(1));
        for (i, &c) in s.iter().enumerate() {
            prefix.push(prefix[i] * base + Mersenne::new(c.into()));
            power.push(power[i] * base);
        }
        Self {
            prefix: prefix.into(),
            power: power.into(),
            base,
        }
    }

    /// The hash `H(s[l..r))` of `range = [l, r)`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `l > r` or `r > n`.
    pub fn hash(&self, range: impl std::ops::RangeBounds<usize>) -> SequenceHash<K> {
        let (l, r) = to_half_open(self.len(), range);
        SequenceHash {
            value: self.prefix[r] - self.prefix[l] * self.power[r - l],
            base: self.base,
            len: r - l,
        }
    }

    /// The length `n` of `s`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.prefix.len() - 1
    }

    /// Whether `n = 0`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        self.prefix.len() == 1
    }
}
