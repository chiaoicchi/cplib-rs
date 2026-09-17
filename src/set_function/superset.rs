use crate::algebra::{Ring, Semiring};

/// The superset zeta transform of `f`, in place.
///
/// # Definition
/// `(Zf)(S) = Σ_{T⊇S} f(T)`, the zeta transform of the boolean lattice on `[n]` ordered by reverse
/// inclusion.
///
/// # Complexity
/// - Time: O(2^n n)
/// - Space: O(1)
///
/// # Panics
/// Panics if `f.len()` is not a power of two.
pub fn superset_zeta<R: Semiring>(ring: &R, f: &mut [R::Value]) {
    assert!(
        f.len().is_power_of_two(),
        "length must be a power of two: len={}",
        f.len()
    );
    let mut w = 1;
    while w < f.len() {
        for block in f.chunks_mut(w << 1) {
            let (lo, hi) = block.split_at_mut(w);
            for (a, b) in lo.iter_mut().zip(hi.iter()) {
                *a = ring.add(a, b);
            }
        }
        w <<= 1;
    }
}

/// The superset Mobius transform of `f`, in place, inverting [`superset_zeta`].
///
/// # Definition
/// `(Z^{-1}f)(S) = Σ_{T⊇S} (-1)^{|T|-|S|} f(T)`, by Mobius inversion on the boolean lattice ordered
/// by reverse inclusion, whose Mobius function is `μ(S, T) = (-1)^{|T|-|S|}`.
///
/// # Complexity
/// - Time: O(2^n n)
/// - Space: O(1)
///
/// # Panics
/// Panics if `f.len()` is not a power of two.
pub fn superset_mobius<R: Ring>(ring: &R, f: &mut [R::Value]) {
    assert!(
        f.len().is_power_of_two(),
        "length must be a power of two: len={}",
        f.len()
    );
    let mut w = 1;
    while w < f.len() {
        for block in f.chunks_mut(w << 1) {
            let (lo, hi) = block.split_at_mut(w);
            for (a, b) in lo.iter_mut().zip(hi.iter()) {
                *a = ring.add(a, &ring.neg(b));
            }
        }
        w <<= 1;
    }
}

/// The intersection convolution of `f` and `g`.
///
/// # Definition
/// `(fg)(S) = Σ_{T∩U=S} f(T) g(U)`, the product of the monoid algebra of the intersection
/// semilattice on `[n]`. Since `T⊇S` and `U⊇S` iff `T∩U⊇S`, [`superset_zeta`] carries it to the
/// pointwise product: `Z(fg) = Z(f) Z(g)`.
///
/// # Complexity
/// - Time: O(2^n n)
/// - Space: O(2^n)
///
/// # Panics
/// Panics if `f` and `g` differ in length.
/// Panics if that length is not a power of two.
pub fn and_convolve<R: Ring>(
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
    superset_zeta(ring, &mut f);
    superset_zeta(ring, &mut g);
    for (fi, gi) in f.iter_mut().zip(g.iter()) {
        *fi = ring.mul(fi, gi);
    }
    superset_mobius(ring, &mut f);
    f
}
