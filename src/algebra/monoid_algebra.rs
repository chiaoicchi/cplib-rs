use crate::algebra::{Field, Ring, RootOfUnity, Semiring};
use crate::arithmetic::dirichlet::dirichlet_convolve;
use crate::arithmetic::gcd::gcd_convolve;
use crate::arithmetic::lcm::lcm_convolve;
use crate::bitwise::and::and_convolve;
use crate::bitwise::or::or_convolve;
use crate::bitwise::xor::xor_convolve;
use crate::cyclic::cyclic_convolve;

/// A convolution over `R`, the type of the products of the named monoid algebras.
pub type Convolve<R> =
    fn(&R, Vec<<R as Semiring>::Value>, Vec<<R as Semiring>::Value>) -> Vec<<R as Semiring>::Value>;

/// The monoid algebra `R[M]` of a monoid `M` of order `order`, as a semiring.
///
/// # Definition
/// `R[M]` is the free `R`-module on `M`, with the product extended bilinearly from `M`:
/// `e_x e_y = e_{xy}`. Its identity is `e_id` for the identity `id` of `M`. Elements are stored as
/// vectors of length `order` with the coefficient of `x` at index `x`, and the product is the
/// convolution `mul`.
///
/// # Contract
/// `mul(ring, f, g)` is the product of `R[M]` on vectors of length `order`, for a monoid `M` on
/// `[0, order)` with identity `id`.
pub struct MonoidAlgebra<R, F> {
    ring: R,
    mul: F,
    id: usize,
    order: usize,
}

impl<R: Semiring, F: Fn(&R, Vec<R::Value>, Vec<R::Value>) -> Vec<R::Value>> MonoidAlgebra<R, F> {
    /// `R[M]` for the monoid `M` on `[0, order)` with product `mul` and identity `id`.
    ///
    /// # Panics
    /// Panics if `id >= order`.
    pub fn new(ring: R, mul: F, id: usize, order: usize) -> Self {
        assert!(
            id < order,
            "id must be less than order: id={id}, order={order}"
        );
        Self {
            ring,
            mul,
            id,
            order,
        }
    }

    /// The order of `M`.
    pub fn order(&self) -> usize {
        self.order
    }
}

impl<R: Ring> MonoidAlgebra<R, Convolve<R>> {
    /// `R[((Z/2)^k, xor)]` on `[0, 2^k)`, whose identity is `0`.
    pub fn xor(ring: R, k: u32) -> Self
    where
        R: Field,
    {
        Self::new(ring, xor_convolve, 0, 1 << k)
    }

    /// `R[(2^[k], ∪)]` on `[0, 2^k)`, whose identity is `0`.
    pub fn or(ring: R, k: u32) -> Self {
        Self::new(ring, or_convolve, 0, 1 << k)
    }

    /// `R[(2^[k], ∩)]` on `[0, 2^k)`, whose identity is `[k]`.
    pub fn and(ring: R, k: u32) -> Self {
        Self::new(ring, and_convolve, (1 << k) - 1, 1 << k)
    }

    /// `R[(N, gcd)]` on `[0, n)`, whose identity is `0`.
    pub fn gcd(ring: R, n: usize) -> Self {
        Self::new(ring, gcd_convolve, 0, n)
    }

    /// `R[(N, lcm)/I]` on `[0, n)`, `I = {0} ∪ [n, ∞)`, whose identity is `1`.
    ///
    /// # Panics
    /// Panics if `n < 2`.
    pub fn lcm(ring: R, n: usize) -> Self {
        Self::new(ring, lcm_convolve, 1, n)
    }
}

impl<R: Semiring> MonoidAlgebra<R, Convolve<R>> {
    /// `R[(N, x)/I]` on `[0, n)`, `I = {0} ∪ [n, ∞)`, whose identity is `1`.
    ///
    /// # Panics
    /// Panics if `n < 2`.
    pub fn dirichlet(ring: R, n: usize) -> Self {
        Self::new(ring, dirichlet_convolve, 1, n)
    }
}

impl<R: RootOfUnity> MonoidAlgebra<R, Convolve<R>> {
    /// `R[Z/nZ]` on `[0, n)`, whose identity is `0`.
    ///
    /// # Panics
    /// Panics if `n` is not a power of two.
    pub fn cyclic(ring: R, n: usize) -> Self {
        Self::new(ring, cyclic_convolve, 0, n)
    }
}

impl<R: Semiring<Value: Clone>, F: Fn(&R, Vec<R::Value>, Vec<R::Value>) -> Vec<R::Value>> Semiring
    for MonoidAlgebra<R, F>
{
    type Value = Vec<R::Value>;

    /// # Complexity
    /// - Time: O(order)
    /// - Space: O(order)
    fn zero(&self) -> Vec<R::Value> {
        (0..self.order).map(|_| self.ring.zero()).collect()
    }

    /// # Complexity
    /// - Time: O(order)
    /// - Space: O(order)
    fn one(&self) -> Vec<R::Value> {
        let mut e = self.zero();
        e[self.id] = self.ring.one();
        e
    }

    /// # Complexity
    /// - Time: O(order)
    /// - Space: O(order)
    ///
    /// # Panics
    /// Panics if the lengths of `a`, `b` and `order` differ.
    fn add(&self, a: &Vec<R::Value>, b: &Vec<R::Value>) -> Vec<R::Value> {
        assert!(
            a.len() == self.order && b.len() == self.order,
            "length mismatch: lhs={}, rhs={}, order={}",
            a.len(),
            b.len(),
            self.order
        );
        a.iter().zip(b).map(|(x, y)| self.ring.add(x, y)).collect()
    }

    /// # Complexity
    /// - Time: that of `mul`, plus O(order)
    /// - Space: that of `mul`, plus O(order)
    ///
    /// # Panics
    /// Panics if the lengths of `a`, `b` and `order` differ.
    fn mul(&self, a: &Vec<R::Value>, b: &Vec<R::Value>) -> Vec<R::Value> {
        assert!(
            a.len() == self.order && b.len() == self.order,
            "length mismatch: lhs={}, rhs={}, order={}",
            a.len(),
            b.len(),
            self.order
        );
        (self.mul)(&self.ring, a.clone(), b.clone())
    }
}

impl<R: Ring<Value: Clone>, F: Fn(&R, Vec<R::Value>, Vec<R::Value>) -> Vec<R::Value>> Ring
    for MonoidAlgebra<R, F>
{
    /// # Complexity
    /// - Time: O(order)
    /// - Space: O(order)
    ///
    /// # Panics
    /// Panics if the length of `a` differ from `order`.
    fn neg(&self, a: &Vec<R::Value>) -> Vec<R::Value> {
        assert!(
            a.len() == self.order,
            "length mismatch: a={}, order={}",
            a.len(),
            self.order
        );
        a.iter().map(|x| self.ring.neg(x)).collect()
    }
}
