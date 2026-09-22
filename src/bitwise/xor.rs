use crate::algebra::{Commutative, Field, Group, Monoid, Ring, Semigroup, Zero};

/// The monoid of `T` under the bitwise xor.
///
/// # Definition
/// Identifying `T` with `{0, 1}^w` by its binary digits, `w` the number of bits of `T`, it is the
/// `w`-fold direct product of `F_2 = ({0, 1}, 0, +)`: `op(a, b) = a ^ b`, and `id()` is `0`.
pub struct Xor<T>(std::marker::PhantomData<T>);
impl<T> Xor<T> {
    pub const fn new() -> Self {
        Self(std::marker::PhantomData)
    }
}
impl<T> Default for Xor<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Clone for Xor<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for Xor<T> {}

impl<T: Clone + std::ops::BitXor<Output = T>> Semigroup for Xor<T> {
    type Value = T;
    fn op(&self, a: &T, b: &T) -> T {
        a.clone() ^ b.clone()
    }
}
impl<T: Clone + std::ops::BitXor<Output = T> + Zero> Monoid for Xor<T> {
    fn id(&self) -> T {
        T::zero()
    }
}
impl<T: Clone + std::ops::BitXor<Output = T> + Zero> Group for Xor<T> {
    fn inv(&self, a: &Self::Value) -> Self::Value {
        a.clone()
    }
}
impl<T: Clone + std::ops::BitOr<Output = T>> Commutative for Xor<T> {}

/// The Walsh-Hadamard transform of `f`, in place.
///
/// # Definition
/// `(Hf)(S) = Σ_T (-1)^{|S∩T|} f(T)`, the expansion of `f` in the characters
/// `T -> (-1)^{|S∩T|}` of the group of the subsets of `[n]` under symmetric difference.
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
/// `(H^{-1}f)(S) = 2^{-n} Σ_T (-1)^{|S∩T|} f(T)`.
///
/// # Contract
/// `2^n` is invertible in `R`, that is the characteristic is not `2`, unless `n = 0`.
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

/// The xor convolution of `f` and `g`.
///
/// # Definition
/// `(fg)(S) = Σ_{T△U=S} f(T) g(U)`, the product of the group algebra of the subsets of `[n]` under
/// symmetric difference.
///
/// # Contract
/// `2^n` is invertible in `R`, that is the characteristic is not `2`, unless `n = 0`.
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
