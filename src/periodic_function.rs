//! Periodic functions on `Z` and the group algebra they form.
//!
//! A periodic function `f` of period `n` is a function on `Z/nZ` and is stored as a slice of length
//! `n` with `f(x)` at index `x`. Throughout this module `n` is that length, a power of two, `ω` is
//! a primitive `n`-th root of unity of `R`, `ω_{2h} = ω^{n/2h}` is the primitive `2h`-th root of
//! unity it determines, and `s`, `x` denote residues modulo `n`.
use crate::algebra::RootOfUnity;

/// The discrete Fourier transform on `f`, in place.
///
/// # Definition
/// `(Ff)(s) = Σ_x ω^{sx} f(x)`, the expansion of `f` in the characters `x -> ω^{sx}` of `Z/nZ`,
/// which exist iff `R` has a primitive `n`-th root of unity.
///
/// # Complexity
/// - Time: O(n log n)
/// - Space: O(n)
///
/// # Panics
/// Panics if `n` is not a power of two.
/// Panics if `R` has no primitive `n`-th root of unity.
pub fn dft<R: RootOfUnity>(ring: &R, f: &mut [R::Value]) {
    let table = twiddles(ring, root_of_unity(ring, f.len()), f.len());
    dif(ring, f, &table);
    bit_reverse(f);
}

/// The inverse discrete Fourier transform of `f`, in place, inverting [`dft`].
///
/// # Definition
/// `(F^{-1}f)(x) = n^{-1} Σ_s ω^{-sx} f(s)`, by the orthogonality of the characters of `Z/nZ`:
/// `Σ_s ω^{s(x-y)}` is `n` if `x = y` and `0` otherwise. A primitive `n`-th root of unity exists
/// only if `n != 0` in `R`, so `n^{-1}` exists whenever `F` does.
///
/// # Complexity
/// - Time: O(n log n)
/// - Space: O(n)
///
/// # Panics
/// Panics if `n` is not a power of two.
/// Panics if `R` has no primitive `n`-th root of unity.
pub fn inverse_dft<R: RootOfUnity>(ring: &R, f: &mut [R::Value]) {
    let omega = root_of_unity(ring, f.len());
    let table = twiddles(ring, ring.inv(&omega), f.len());
    bit_reverse(f);
    dit(ring, f, &table);
    let mut n = ring.one();
    for _ in 0..f.len().trailing_zeros() {
        n = ring.add(&n, &n);
    }
    let inv_n = ring.inv(&n);
    for x in f.iter_mut() {
        *x = ring.mul(x, &inv_n);
    }
}

/// The cyclic convolution of `f` and `g`.
///
/// # Definition
/// `(fg)(x) = Σ_{y+z=x} f(y) g(z)`, the product of the group algebra `R[Z/nZ] = R[t]/(t^n - 1)`.
/// The characters are multiplicative, `ω^{s(y+z)} = ω^{sy} ω^{sz}`, so [`dft`] carries it to the
/// pointwise product: `F(fg) = F(f) F(g)`.
///
/// # Complexity
/// - Time: O(n log n)
/// - Space: O(n)
///
/// # Panics
/// Panics if `f` and `g` differ in length.
/// Panics if that length is not a power of two.
/// Panics if `R` has no primitive `n`-th root of unity.
pub fn cyclic_convolve<R: RootOfUnity>(
    ring: &R,
    mut f: Vec<R::Value>,
    mut g: Vec<R::Value>,
) -> Vec<R::Value> {
    assert_eq!(
        f.len(),
        g.len(),
        "f and g must have the same length: f={}, g={}",
        f.len(),
        g.len(),
    );
    let n = f.len();
    let omega = root_of_unity(ring, n);
    let inv_table = twiddles(ring, ring.inv(&omega), n);
    let table = twiddles(ring, omega, n);
    dif(ring, &mut f, &table);
    dif(ring, &mut g, &table);
    for (fi, gi) in f.iter_mut().zip(g.iter()) {
        *fi = ring.mul(fi, gi);
    }
    dit(ring, &mut f, &inv_table);
    let mut n = ring.one();
    for _ in 0..f.len().trailing_zeros() {
        n = ring.add(&n, &n);
    }
    let inv_n = ring.inv(&n);
    for x in f.iter_mut() {
        *x = ring.mul(x, &inv_n);
    }
    f
}

/// A primitive `n`-th root of unity of `R`.
///
/// # Panics
/// Panics if `R` has no primitive `n`-th root of unity.
fn root_of_unity<R: RootOfUnity>(ring: &R, n: usize) -> R::Value {
    ring.root_of_unity(n)
        .unwrap_or_else(|| panic!("no primitive {n}-th root of unity"))
}

/// The powers of `omega` at every stage.
///
/// # Definition
/// `ω_{2h}^j` at index `h + j` for `j < h` and `h = n/2, ..., 1`, where `ω` is `omega`.
///
/// # Contract
/// `omega` is a primitive `n`-th root of unity.
///
/// # Complexity
/// - Time: O(n)
/// - Space: O(n)
///
/// # Panics
/// Panics if `n` is not a power of two.
pub(crate) fn twiddles<R: RootOfUnity>(ring: &R, mut omega: R::Value, n: usize) -> Vec<R::Value> {
    assert!(n.is_power_of_two(), "n must be a power of two: n={n}");
    let mut table: Vec<R::Value> = (0..n).map(|_| ring.one()).collect();
    let mut h = n >> 1;
    while h >= 1 {
        for j in 1..h {
            table[h + j] = ring.mul(&table[h + j - 1], &omega);
        }
        omega = ring.mul(&omega, &omega);
        h >>= 1;
    }
    table
}

/// The stages of [`dft`] from natural order to bit-reversed order, in place.
///
/// # Definition
/// For `h = n/2, ..., 1`, `(a, b) -> (a + b, (a - b) ω_{2h}^j)` on the pairs `(f(x), f(x + h))`,
/// `j = x mod h`, of each block of length `2h`: the reduction of the block modulo `t^h - 1` and
/// modulo `t^h + 1`, the latter followed by `t -> ω_{2h} u`, since
/// `R[t]/(t^{2h} - 1) = R[t]/(t^h - 1) x R[t]/(t^h + 1)`.
///
/// # Contract
/// `table` is [`twiddles`] of `n` for the root of unity of the transform.
///
/// # Complexity
/// - Time: O(n log n)
/// - Space: O(1)
pub(crate) fn dif<R: RootOfUnity>(ring: &R, f: &mut [R::Value], table: &[R::Value]) {
    let mut h = f.len() >> 1;
    while h >= 1 {
        for block in f.chunks_mut(h << 1) {
            let (lo, hi) = block.split_at_mut(h);
            for ((a, b), omega) in lo.iter_mut().zip(hi.iter_mut()).zip(&table[h..h << 1]) {
                (*a, *b) = (ring.add(a, b), ring.mul(&ring.add(a, &ring.neg(b)), omega));
            }
        }
        h >>= 1;
    }
}

/// The stages of [`inverse_dft`] from bit-reversed order to natural order, in place.
///
/// # Definition
/// For `h = 1, ..., n/2`, `(a, b) -> (a + ω_{2h}^j b, a - ω_{2h}^j b)` on the pairs
/// `(f(x), f(x + h))`, `j = x mod h`, of each block of length `2h`: the transform of a block from
/// the transforms of its even and odd parts, since `Σ_y ω_{2h}^{sy} f(y)` splits by the parity of
/// `y` into `Σ_z ω_h^{sz} f(2z) + ω_{2h}^s Σ_z ω_h^{sz} f(2z + 1)`.
///
/// # Contract
/// `table` is [`twiddles`] of `n` for the root of unity of the transform.
///
/// # Complexity
/// - Time: O(n log n)
/// - Space: O(1)
pub(crate) fn dit<R: RootOfUnity>(ring: &R, f: &mut [R::Value], table: &[R::Value]) {
    let mut h = 1;
    while h < f.len() {
        for block in f.chunks_mut(h << 1) {
            let (lo, hi) = block.split_at_mut(h);
            for ((a, b), omega) in lo.iter_mut().zip(hi.iter_mut()).zip(&table[h..h << 1]) {
                let t = ring.mul(b, omega);
                (*a, *b) = (ring.add(a, &t), ring.add(a, &ring.neg(&t)));
            }
        }
        h <<= 1;
    }
}

/// Permutes `f` by the bit reversal of the index, in place.
///
/// # Definition
/// `f(x) -> f(rev(x))`, where `rev` reverses the `log n` bits of `x`. It is an involution.
///
/// # Complexity
/// - Time: O(n)
/// - Space: O(1)
fn bit_reverse<T>(f: &mut [T]) {
    let n = f.len();
    let k = n.trailing_zeros();
    for i in 1..n {
        let j = i.reverse_bits() >> (usize::BITS - k);
        if i < j {
            f.swap(i, j);
        }
    }
}
