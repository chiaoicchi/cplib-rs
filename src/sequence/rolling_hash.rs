use crate::algebra::{Monoid, Semigroup};
use crate::num::mersenne::Mersenne;
use crate::range::to_half_open;

/// A polynomial hash of a sequence together with its length.
///
/// # Definition
/// For a base `b` in `Z/MZ` with `M = 2^K - 1`, the hash of a sequence `t` of length `m` is
/// `H(t) = t[0] b^(m-1) + t[1] b^(m-2) + ... + t[m-1]` in `Z/MZ`. Two hashes are equal iff
/// their lengths, values and base are equal. For `t != t'` of length at most `n`, `H(t) = H(t')`
/// with probability at most `n / M` over a uniformly random `b`, since `H(t) - H(t')` is a
/// nonzero polynomial in `base` of degree less than `n`.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct SequenceHash<const K: u32> {
    value: Mersenne<K>,
    base: Mersenne<K>,
    len: usize,
}

impl<const K: u32> SequenceHash<K> {
    /// Returns the hash value.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn value(self) -> Mersenne<K> {
        self.value
    }

    /// Returns the length of the hashed sequence.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(self) -> usize {
        self.len
    }

    /// Returns `true` if the hashed sequence is empty.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(self) -> bool {
        self.len == 0
    }
}

/// The monoid of sequence hashes under concatenation.
///
/// # Definition
/// `H(t u) = H(t) b^|u| + H(u)`, with the empty sequence as the identity.
#[derive(Clone, Copy)]
pub struct Concat<const K: u32>(pub Mersenne<K>);
impl<const K: u32> Semigroup for Concat<K> {
    type Value = SequenceHash<K>;
    /// # Complexity
    /// - Time: O(log b.len())
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `a.base` and `b.base` are not same.
    fn op(&self, a: &SequenceHash<K>, b: &SequenceHash<K>) -> SequenceHash<K> {
        assert!(
            self.0 == a.base,
            "bases differ: expected={}, a={}",
            self.0,
            a.base
        );
        assert!(
            self.0 == b.base,
            "bases differ: expected={}, b={}",
            self.0,
            b.base
        );
        SequenceHash {
            value: a.value * self.0.pow(b.len() as u64) + b.value,
            base: self.0,
            len: a.len + b.len,
        }
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

/// The prefix hashes of a sequence, giving the hash of any contiguous subsequence in O(1).
///
/// # Definition
/// For a sequence `s` of length `n` and a base `b`, `prefix[i] = H(s[0..i))` for `i` in `[0, n]`
/// and `power[i] = b^i`. Then `H(s[l..r)) = prefix[r] - prefix[l] b^(r-l)`.
///
/// # Contract
/// - `base` is chosen uniformly at random from `Z/MZ` at run time; the collision bound of
///   `SequenceHash` holds only for a random base, and a fixed base admits adversarial inputs.
/// - Hashes are comparable only if they were computed with the same base.
///
/// # Invariants
/// `prefix.len() = power.len() = n + 1`.
///
/// # Complexity
/// - Space: O(n)
pub struct RollingHash<const K: u32> {
    prefix: Box<[Mersenne<K>]>,
    power: Box<[Mersenne<K>]>,
    base: Mersenne<K>,
}

impl<const K: u32> RollingHash<K> {
    /// Constructs the prefix hashes of `s` with base `base`.
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

    /// Returns the hash of `s[range]`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `range` is out of bounds or `l > r`.
    pub fn hash(&self, range: impl std::ops::RangeBounds<usize>) -> SequenceHash<K> {
        let (l, r) = to_half_open(self.len(), range);
        assert!(
            l <= r,
            "left bound must be less than or equal to right bound: l={l}, r={r}"
        );
        assert!(r <= self.len(), "range out of bounds: range=[{l}, {r})");
        SequenceHash {
            value: self.prefix[r] - self.prefix[l] * self.power[r - l],
            base: self.base,
            len: r - l,
        }
    }

    /// Returns the length of `s`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> usize {
        self.prefix.len() - 1
    }

    /// Returns `true` if `s` is empty.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        self.prefix.len() == 1
    }
}
