use crate::algebra::{Commutative, Idempotent, Monoid, Ring, Semigroup, Semiring, Zero};

/// The monoid of `T` under the bitwise or.
///
/// # Definition
/// Identifying `T` with `{0, 1}^w` by its binary digits, `w` the number of bits of `T`, it is the
/// `w`-fold direct product of `({0, 1}, 0, v)`: `op(a, b) = a | b`, and `id()` is `0`.
pub struct Or<T>(std::marker::PhantomData<T>);
impl<T> Or<T> {
    pub const fn new() -> Self {
        Self(std::marker::PhantomData)
    }
}
impl<T> Default for Or<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Clone for Or<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for Or<T> {}

impl<T: Clone + std::ops::BitOr<Output = T>> Semigroup for Or<T> {
    type Value = T;
    fn op(&self, a: &T, b: &T) -> T {
        a.clone() | b.clone()
    }
}
impl<T: Clone + std::ops::BitOr<Output = T> + Zero> Monoid for Or<T> {
    fn id(&self) -> T {
        T::zero()
    }
}
impl<T: Clone + std::ops::BitOr<Output = T>> Commutative for Or<T> {}
impl<T: Clone + std::ops::BitOr<Output = T>> Idempotent for Or<T> {}

/// The subset zeta transform of `f`, in place.
///
/// # Definition
/// `(Zf)(S) = Σ_{T⊆S} f(T)`, the zeta trnsform of the boolean lattice on `[n]` ordered by
/// inclusion.
///
/// # Complexity
/// - Time: O(2^n n)
/// - Space: O(1)
///
/// # Panics
/// Panics if `f.len()` is not a power of two.
pub fn subset_zeta<R: Semiring>(ring: &R, f: &mut [R::Value]) {
    assert!(
        f.len().is_power_of_two(),
        "length must be a power of two: len={}",
        f.len()
    );
    let mut w = 1;
    while w < f.len() {
        for block in f.chunks_mut(w << 1) {
            let (lo, hi) = block.split_at_mut(w);
            for (a, b) in hi.iter_mut().zip(lo.iter()) {
                *a = ring.add(a, b);
            }
        }
        w <<= 1;
    }
}

/// The subset Mobius transform of `f`, in place, inverting [`subset_zeta`].
///
/// # Definition
/// `(Z^{-1}f)(S) = Σ_{T⊆S} (-1)^{|S|-|T|} f(T)`.
///
/// # Complexity
/// - Time: O(2^n n)
/// - Space: O(1)
///
/// # Panics
/// Panics if `f.len()` is not a power of two.
pub fn subset_mobius<R: Ring>(ring: &R, f: &mut [R::Value]) {
    assert!(
        f.len().is_power_of_two(),
        "length must be a power of two: len={}",
        f.len()
    );
    let mut w = 1;
    while w < f.len() {
        for block in f.chunks_mut(w << 1) {
            let (lo, hi) = block.split_at_mut(w);
            for (a, b) in hi.iter_mut().zip(lo.iter()) {
                *a = ring.add(a, &ring.neg(b));
            }
        }
        w <<= 1;
    }
}

/// The or convolution of `f` and `g`.
///
/// # Definition
/// `(fg)(S) = Σ_{T∪U=S} f(T) g(U)`, the product of the monoid algebra of the subsets of `[n]` under
/// union.
///
/// # Complexity
/// - Time: O(2^n n)
/// - Space: O(2^n)
///
/// # Panics
/// Panics if `f` and `g` differ in length.
/// Panics if that length is not a power of two.
pub fn or_convolve<R: Ring>(ring: &R, mut f: Vec<R::Value>, mut g: Vec<R::Value>) -> Vec<R::Value> {
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
    subset_zeta(ring, &mut f);
    subset_zeta(ring, &mut g);
    for (fi, gi) in f.iter_mut().zip(g.iter()) {
        *fi = ring.mul(fi, gi);
    }
    subset_mobius(ring, &mut f);
    f
}
