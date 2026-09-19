use crate::algebra::{Field, Ring};

/// The Walsh-Hadamard transform of `f`, in place.
///
/// # Definition
/// `(Hf)(S) = Σ_T (-1)^{|S∩T|} f(T)`, the expansion of `f` in the characters of the group
/// `(Z/2)^n`, where the character indexed by `S` sends `T` to `(-1)^{|S∩T|}`. It is its own inverse
/// up to a scalar: `H(Hf) = 2^n f`.
///
/// # Complexity
/// - Time: O(2^n n)
/// - Space: O(1)
///
/// # Panics
/// Panics if `f.len()` is not a power of two.
pub fn walsh_hadamard<R: Ring>(ring: &R, f: &mut [R::Value]) {
    assert!(
        f.len().is_power_of_two(),
        "length must be a power of two: len={}",
        f.len()
    );
    let mut w = 1;
    while w < f.len() {
        for block in f.chunks_mut(w << 1) {
            let (lo, hi) = block.split_at_mut(w);
            for (a, b) in lo.iter_mut().zip(hi.iter_mut()) {
                (*a, *b) = (ring.add(a, b), ring.add(a, &ring.neg(b)));
            }
        }
        w <<= 1;
    }
}

/// The inverse Walsh-Hadamard transform of `f`, in place, inverting [`walsh_hadamard`].
///
/// # Definition
/// `(H^{-1}f)(S) = 2^{-n} Σ_T (-1)^{|S∩T|} f(T)`, by the orthogonality of the characters of
/// `(Z/2)^n`: `Σ_T (-1)^{|S∩T|} (-1)^{|U∩T|}` is `2^n` if `S = U` and `0` otherwise.
///
/// # Contract
/// `2^n` is invertible in `R`, that is the characteristic is `0` or odd.
///
/// # Complexity
/// - Time: O(2^n n)
/// - Space: O(1)
///
/// # Panics
/// Panics if `f.len()` is not a power of two.
pub fn inverse_walsh_hadamard<R: Field>(ring: &R, f: &mut [R::Value]) {
    walsh_hadamard(ring, f);
    let mut c = ring.one();
    let two = ring.add(&ring.one(), &ring.one());
    for _ in 0..f.len().trailing_zeros() {
        c = ring.mul(&c, &two);
    }
    let c = ring.inv(&c);
    for x in f.iter_mut() {
        *x = ring.mul(&c, x);
    }
}

/// The symmetric difference convolution of `f` and `g`.
///
/// # Definition
/// `(fg)(S) = Σ_{T△U=S} f(T) g(U)`, the product of the group algebra of `(Z/2)^n`. The characters
/// are multiplicative, `(-1)^{|S∩(T△U)|} = (-1)^{|S∩T|} (-1)^{|S∩U|}`, so [`walsh_hadamard`]
/// carries the product to the pointwise product: `H(fg) = H(f) H(g)`.
///
/// # Contract
/// `2^n` is invertible in `R`, that is the characteristic is `0` or odd.
///
/// # Complexity
/// - Time: O(2^n n)
/// - Space: O(2^n)
///
/// # Panics
/// Panics if `f` and `g` differ in length.
/// Panics if that length is not a power of two.
pub fn xor_convolve<R: Field>(
    ring: &R,
    mut f: Vec<R::Value>,
    mut g: Vec<R::Value>,
) -> Vec<R::Value> {
    assert_eq!(
        f.len(),
        g.len(),
        "f and g must have the same length: f={}, g={}",
        f.len(),
        g.len()
    );
    assert!(
        f.len().is_power_of_two(),
        "length must be a power of two: len={}",
        f.len()
    );
    walsh_hadamard(ring, &mut f);
    walsh_hadamard(ring, &mut g);
    for (fi, gi) in f.iter_mut().zip(g.iter()) {
        *fi = ring.mul(fi, gi);
    }
    inverse_walsh_hadamard(ring, &mut f);
    f
}
