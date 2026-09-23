/// A multiset `S` of `[0, 2^bits)`, as a binary trie.
///
/// # Definition
/// `n = |S|` counted with multiplicity, and `bits` is the number of bits of its elements. `flip` is
/// the value xored into every element by [`xor_all`](Self::xor_all), so that the element stored as
/// `y` is `y ^ flip`.
///
/// # Invariants
/// - `child[v]` are the two children of the node `v`, or `u32::MAX` if absent, and `count[v]` is
///   the number of stored values passing through `v`, with the root at `0`.
/// - Every node reachable from the root has `count[v] > 0`, `count[0] = n`, and `distinct` is the
///   number of distinct elements.
///
/// # Complexity
/// - Space: O(k bits), where `k` is the number of distinct elements
pub struct BinaryTrie {
    child: Vec<[u32; 2]>,
    count: Vec<u64>,
    distinct: u64,
    flip: u64,
    bits: u32,
}

impl BinaryTrie {
    /// The empty multiset of `[0, 2^bits)`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `bits > 64`.
    pub fn new(bits: u32) -> Self {
        assert!(bits <= 64, "bits must be at most 64: bits={bits}");
        Self {
            child: vec![[!0; 2]],
            count: vec![0],
            distinct: 0,
            flip: 0,
            bits,
        }
    }

    /// The multiset of `[0, 2^bits)` of the elements of `xs`.
    ///
    /// # Complexity
    /// - Time: O(|xs| log |xs| + k bits), where `k` is the number of distinct element
    /// - Space: O(k bits)
    ///
    /// # Panics
    /// Panics if `bits > 64`, or `xs[i] >= 2^bits` for some `i`.
    pub fn from_vec(bits: u32, mut xs: Vec<u64>) -> Self {
        let mut trie = Self::new(bits);
        xs.sort_unstable();
        let mut parent = vec![!0u32];
        let mut path = vec![0u32; bits as usize + 1];
        let mut previous = None;
        for x in xs {
            assert!(
                bits == 64 || x < 1 << bits,
                "out of bounds: x={x}, bits={bits}",
            );
            if previous == Some(x) {
                trie.count[path[bits as usize] as usize] += 1;
                continue;
            }
            trie.distinct += 1;
            let start = match previous {
                None => 0,
                Some(p) => bits - (64 - (p ^ x).leading_zeros()),
            };
            let mut v = path[start as usize];
            for h in (0..bits - start).rev() {
                let b = (x >> h & 1) as usize;
                let c = trie.child.len() as u32;
                trie.child[v as usize][b] = c;
                trie.child.push([!0; 2]);
                trie.count.push(0);
                parent.push(v);
                v = c;
                path[(bits - h) as usize] = c;
            }
            trie.count[v as usize] = 1;
            previous = Some(x);
        }
        for v in (1..trie.count.len()).rev() {
            trie.count[parent[v] as usize] += trie.count[v];
        }
        trie
    }

    /// Inserts `c` copies of `x` into `S`, and returns the number of copies of `x` afterwards.
    ///
    /// # Complexity
    /// - Time: O(bits)
    /// - Space: O(bits)
    ///
    /// # Panics
    /// Panics if `x >= 2^bits`.
    pub fn insert_count(&mut self, x: u64, c: u64) -> u64 {
        assert!(
            self.bits == 64 || x < 1 << self.bits,
            "out of bounds: x={x}, bits={}",
            self.bits,
        );
        if c == 0 {
            return self.count(x);
        }
        let y = x ^ self.flip;
        let mut v = 0;
        self.count[v] += c;
        for h in (0..self.bits).rev() {
            let b = (y >> h & 1) as usize;
            if self.child[v][b] == !0 {
                self.child[v][b] = self.child.len() as u32;
                self.child.push([!0; 2]);
                self.count.push(0);
            }
            v = self.child[v][b] as usize;
            self.count[v] += c;
        }
        if self.count[v] == c {
            self.distinct += 1;
        }
        self.count[v]
    }

    /// Inserts `x` into `S`, and returns the number of copies of `x` afterwards.
    ///
    /// # Complexity
    /// - Time: O(bits)
    /// - Space: O(bits)
    ///
    /// # Panics
    /// Panics if `x >= 2^bits`.
    pub fn insert(&mut self, x: u64) -> u64 {
        self.insert_count(x, 1)
    }

    /// Removes at most `c` copies of `x` from `S`, and returns the number of copies removed.
    ///
    /// # Complexity
    /// - Time: O(bits)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x >= 2^bits`.
    pub fn remove_count(&mut self, x: u64, c: u64) -> u64 {
        let c = self.count(x).min(c);
        if c == 0 {
            return 0;
        }
        if self.count(x) == c {
            self.distinct -= 1;
        }
        let y = x ^ self.flip;
        let mut v = 0;
        self.count[v] -= c;
        for h in (0..self.bits).rev() {
            let b = (y >> h & 1) as usize;
            let u = self.child[v][b] as usize;
            self.count[u] -= c;
            if self.count[u] == 0 {
                self.child[v][b] = !0;
                break;
            }
            v = u;
        }
        c
    }

    /// Removes one copy of `x` from `S`, and returns the number of copies removed.
    ///
    /// # Complexity
    /// - Time: O(bits)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x >= 2^bits`.
    pub fn remove(&mut self, x: u64) -> u64 {
        self.remove_count(x, 1)
    }

    /// Removes every copy of `x` from `S`, and returns the number of copies removed.
    ///
    /// # Complexity
    /// - Time: O(bits)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x >= 2^bits`.
    pub fn remove_all(&mut self, x: u64) -> u64 {
        self.remove_count(x, !0)
    }

    /// The number of copies of `x` in `S`.
    ///
    /// # Complexity
    /// - Time: O(bits)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x >= 2^bits`.
    pub fn count(&self, x: u64) -> u64 {
        assert!(
            self.bits == 64 || x < 1 << self.bits,
            "out of bounds: x={x}, bits={}",
            self.bits
        );
        let y = x ^ self.flip;
        let mut v = 0;
        for h in (0..self.bits).rev() {
            let c = self.child[v][(y >> h & 1) as usize];
            if c == !0 {
                return 0;
            }
            v = c as usize;
        }
        self.count[v]
    }

    /// Whether `x` is in `S`.
    ///
    /// # Complexity
    /// - Time: O(bits)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x >= 2^bits`.
    pub fn contains(&self, x: u64) -> bool {
        self.count(x) > 0
    }

    /// Replaces every element `x` of `S` by `x ^ v`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `v >= 2^bits`.
    pub fn xor_all(&mut self, v: u64) {
        assert!(
            self.bits == 64 || v < 1 << self.bits,
            "out of bounds: v={v}, bits={}",
            self.bits,
        );
        self.flip ^= v;
    }

    /// The `k`-th least element of `{x ^ v: x ∈ S}` where `k` is `0`-indexed, or
    /// `None` if `k >= n`.
    ///
    /// # Complexity
    /// - Time: O(bits)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `v >= 2^bits`.
    pub fn kth_xor(&self, mut k: u64, v: u64) -> Option<u64> {
        assert!(
            self.bits == 64 || v < 1 << self.bits,
            "out of bounds: v={v}, bits={}",
            self.bits,
        );
        if k >= self.len() {
            return None;
        }
        let f = self.flip ^ v;
        let mut node = 0;
        let mut x = 0;
        for h in (0..self.bits).rev() {
            let b = (f >> h & 1) as usize;
            let c = self.child[node][b];
            let cnt = if c == !0 { 0 } else { self.count[c as usize] };
            let b = if k < cnt {
                b
            } else {
                k -= cnt;
                b ^ 1
            };
            node = self.child[node][b] as usize;
            x = x << 1 | (b ^ (f >> h & 1) as usize) as u64;
        }
        Some(x)
    }

    /// The `k`-th least element of `S` where `k` is `0`-indexed, or `None` if `k >= n`.
    ///
    /// # Complexity
    /// - Time: O(bits)
    /// - Space: O(1)
    pub fn kth(&self, k: u64) -> Option<u64> {
        self.kth_xor(k, 0)
    }

    /// The least element of `{x ^ v: x ∈ S}`, or `None` if `S` is empty.
    ///
    /// # Complexity
    /// - Time: O(bits)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `v >= 2^bits`.
    pub fn min_xor(&self, v: u64) -> Option<u64> {
        self.kth_xor(0, v)
    }

    /// The greatest element of `{x ^ v: x ∈ S}`, or `None` if `S` is empty.
    ///
    /// # Complexity
    /// - Time: O(bits)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `v >= 2^bits`.
    pub fn max_xor(&self, v: u64) -> Option<u64> {
        self.kth_xor(self.len().checked_sub(1)?, v)
    }

    /// The least element of `S`, or `None` if `S` is empty.
    ///
    /// # Complexity
    /// - Time: O(bits)
    /// - Space: O(1)
    pub fn min(&self) -> Option<u64> {
        self.kth(0)
    }

    /// The greatest element of `S`, or `None` if `S` is empty.
    ///
    /// # Complexity
    /// - Time: O(bits)
    /// - Space: O(1)
    pub fn max(&self) -> Option<u64> {
        self.kth(self.len().checked_sub(1)?)
    }

    /// The number of elements of `S` less than `x`, counted with multiplicity.
    ///
    /// # Complexity
    /// - Time: O(bits)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x > 2^bits`.
    pub fn rank(&self, x: u64) -> u64 {
        assert!(
            self.bits == 64 || x <= 1 << self.bits,
            "out of bounds: x={x}, bits={}",
            self.bits,
        );
        if self.bits < 64 && x == 1 << self.bits {
            return self.len();
        }
        let mut node = 0;
        let mut r = 0;
        for h in (0..self.bits).rev() {
            let f = (self.flip >> h & 1) as usize;
            let b = (x >> h & 1) as usize;
            if b == 1 {
                let c = self.child[node][f];
                if c != !0 {
                    r += self.count[c as usize];
                }
            }
            let c = self.child[node][f ^ b];
            if c == !0 {
                return r;
            }
            node = c as usize;
        }
        r
    }

    /// The number of elements of `S` in `[l, r)`, counted with multiplicity.
    ///
    /// # Complexity
    /// - Time: O(bits)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `l > r` or `r > 2^bits`.
    pub fn range_count(&self, l: u64, r: u64) -> u64 {
        assert!(
            l <= r,
            "left bound must be at most right bound: l={l}, r={r}"
        );
        self.rank(r) - self.rank(l)
    }

    /// The least element of `S` at least `x`, or `None` if there is none.
    ///
    /// # Complexity
    /// - Time: O(bits)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x > 2^bits`.
    pub fn ceil(&self, x: u64) -> Option<u64> {
        assert!(
            self.bits == 64 || x <= 1 << self.bits,
            "out of bounds: x={x}, bits={}",
            self.bits,
        );
        if self.bits < 64 && x == 1 << self.bits {
            return None;
        }
        let mut node = 0;
        let mut prefix = 0;
        let mut candidate = None;
        for h in (0..self.bits).rev() {
            let f = (self.flip >> h & 1) as usize;
            let b = (x >> h & 1) as usize;
            if b == 0 {
                let c = self.child[node][f ^ 1];
                if c != !0 {
                    candidate = Some((c as usize, prefix << 1 | 1, h));
                }
            }
            let c = self.child[node][f ^ b];
            if c == !0 {
                let (node, prefix, h) = candidate?;
                return Some(self.min_in(node, prefix, h));
            }
            node = c as usize;
            prefix = prefix << 1 | b as u64;
        }
        Some(x)
    }

    /// The greatest element of `S` at most `x`, or `None` if there is none.
    ///
    /// # Complexity
    /// - Time: O(bits)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `x >= 2^bits`.
    pub fn floor(&self, x: u64) -> Option<u64> {
        assert!(
            self.bits == 64 || x < 1 << self.bits,
            "out of bounds: x={x}, bits={}",
            self.bits,
        );
        let mut node = 0;
        let mut prefix = 0;
        let mut candidate = None;
        for h in (0..self.bits).rev() {
            let f = (self.flip >> h & 1) as usize;
            let b = (x >> h & 1) as usize;
            if b == 1 {
                let c = self.child[node][f];
                if c != !0 {
                    candidate = Some((c as usize, prefix << 1, h));
                }
            }
            let c = self.child[node][f ^ b];
            if c == !0 {
                let (node, prefix, h) = candidate?;
                return Some(self.max_in(node, prefix, h));
            }
            node = c as usize;
            prefix = prefix << 1 | b as u64;
        }
        Some(x)
    }

    /// The least element of the subtree of the node `v`, whose values have the `bits - h` high bits
    /// `prefix`.
    ///
    /// # Complexity
    /// - Time: O(h)
    /// - Space: O(1)
    fn min_in(&self, mut v: usize, mut prefix: u64, h: u32) -> u64 {
        for g in (0..h).rev() {
            let f = (self.flip >> g & 1) as usize;
            let b = if self.child[v][f] == !0 { 1 } else { 0 };
            v = self.child[v][f ^ b] as usize;
            prefix = prefix << 1 | b as u64;
        }
        prefix
    }

    /// The greatest element of the subtree of the node `v`, whose values have the `bits - h` high
    /// bits `prefix`.
    ///
    /// # Complexity
    /// - Time: O(h)
    /// - Space: O(1)
    fn max_in(&self, mut v: usize, mut prefix: u64, h: u32) -> u64 {
        for g in (0..h).rev() {
            let f = (self.flip >> g & 1) as usize;
            let b = if self.child[v][f ^ 1] == !0 { 0 } else { 1 };
            v = self.child[v][f ^ b] as usize;
            prefix = prefix << 1 | b as u64;
        }
        prefix
    }

    /// The size `n` of `S`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn len(&self) -> u64 {
        self.count[0]
    }

    /// The number `k` of distinct elements of `S`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn distinct_len(&self) -> u64 {
        self.distinct
    }

    /// Whether `n = 0`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The number `bits` of bits.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn bits(&self) -> u32 {
        self.bits
    }
}
