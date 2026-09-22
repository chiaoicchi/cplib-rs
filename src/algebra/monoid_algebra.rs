use crate::algebra::{Field, Ring, RootOfUnity, Semiring};
use crate::arithmetic::dirichlet::dirichlet_convolve;
use crate::arithmetic::gcd::gcd_convolve;
use crate::arithmetic::lcm::lcm_convolve;
use crate::bitwise::and::and_convolve;
use crate::bitwise::or::or_convolve;
use crate::bitwise::set_power_series::sps_convolve;
use crate::bitwise::xor::xor_convolve;
use crate::cyclic::cyclic_convolve;
use crate::fps::fps_convolve;
use crate::poly::poly_convolve;

/// A convolution over `R`, the type of the products of the named algebras.
pub type Convolve<R> =
    fn(&R, Vec<<R as Semiring>::Value>, Vec<<R as Semiring>::Value>) -> Vec<<R as Semiring>::Value>;

/// The monoid algebra `R[M]` of a monoid `M`, or its quotient `R[M]/R[I]` by the span of an ideal
/// `I` of `M`, as a samiring.
///
/// # Definition
/// `R[M]` is the free `R`-module on `M`, with the product `e_x e_y = e_{xy}` extended bilinearly
/// and the identity `e_id`. The elements of `M`, or of `M \ I`, are indexed by `[0, order)` if
/// `order` is `Some`, and by `N` if it is `None`. An element `Σ_x f(x) e_x` is stored as a vector
/// `f(x) = 0` beyond it, if `order` is `None`.
pub struct MonoidAlgebra<R, F> {
    ring: R,
    mul: F,
    id: usize,
    order: Option<usize>,
}
impl<R: Clone, F: Clone> Clone for MonoidAlgebra<R, F> {
    fn clone(&self) -> Self {
        Self {
            ring: self.ring.clone(),
            mul: self.mul.clone(),
            id: self.id,
            order: self.order,
        }
    }
}
impl<R: Copy, F: Copy> Copy for MonoidAlgebra<R, F> {}

impl<R: Semiring, F: Fn(&R, Vec<R::Value>, Vec<R::Value>) -> Vec<R::Value>> MonoidAlgebra<R, F> {
    /// The algebra with the product `mul` and the identity `e_id`.
    ///
    /// # Contract
    /// `mul` is the product of the algebra on stored vectors, and `e_id` is its identity.
    ///
    /// # Panics
    /// Panics if `order` is `Some` and `id >= order`.
    pub fn new(ring: R, mul: F, id: usize, order: Option<usize>) -> Self {
        if let Some(order) = order {
            assert!(
                id < order,
                "id must be less than order: id={id}, order={order}"
            );
        }
        Self {
            ring,
            mul,
            id,
            order,
        }
    }
    /// The number of basis elements, `None` if it is infinite.
    pub fn order(&self) -> Option<usize> {
        self.order
    }
}

impl<R: Ring> MonoidAlgebra<R, Convolve<R>> {
    /// `R[(Z/2)^k]`, the set functions on `[k]` under [`xor_convolve`].
    pub fn xor(ring: R, k: u32) -> Self
    where
        R: Field,
    {
        Self::new(ring, xor_convolve, 0, Some(1 << k))
    }

    /// The set functions on `[k]` under [`or_convolve`].
    pub fn or(ring: R, k: u32) -> Self {
        Self::new(ring, or_convolve, 0, Some(1 << k))
    }

    /// The set functions on `[k]` under [`and_convolve`].
    pub fn and(ring: R, k: u32) -> Self {
        Self::new(ring, and_convolve, (1 << k) - 1, Some(1 << k))
    }

    /// The set power series on `[k]`, under [`sps_convolve`].
    pub fn sps(ring: R, k: u32) -> Self {
        Self::new(ring, sps_convolve, 0, Some(1 << k))
    }

    /// The arithmetic functions on `[0, n)` under [`gcd_convolve`].
    ///
    /// # Panics
    /// Panics if `n = 0`.
    pub fn gcd(ring: R, n: usize) -> Self {
        Self::new(ring, gcd_convolve, 0, Some(n))
    }

    /// The arithmetic functions on `[0, n)` under [`lcm_convolve`].
    ///
    /// # Panics
    /// Panics if `n <= 1`
    pub fn lcm(ring: R, n: usize) -> Self {
        Self::new(ring, lcm_convolve, 1, Some(n))
    }
}
impl<R: Semiring> MonoidAlgebra<R, Convolve<R>> {
    /// The arithmetic functions on `[0, n)` under [`dirichlet_convolve`].
    ///
    /// # Panics
    /// Panics if `n <= 1`.
    pub fn dirichlet(ring: R, n: usize) -> Self {
        Self::new(ring, dirichlet_convolve, 1, Some(n))
    }
}
impl<R: RootOfUnity> MonoidAlgebra<R, Convolve<R>> {
    /// `R[Z/nZ]`, the functions on `Z/nZ` under [`cyclic_convolve`].
    ///
    /// # Panics
    /// Panics if `n = 0`.
    pub fn cyclic(ring: R, n: usize) -> Self {
        Self::new(ring, cyclic_convolve, 0, Some(n))
    }

    /// `R[x]/(x^n)`, the formal power series of precision `n` under [`fps_convolve`].
    ///
    /// # Panics
    /// Panics if `n = 0`.
    pub fn fps(ring: R, n: usize) -> Self {
        Self::new(
            ring,
            |ring, f, g| {
                let n = f.len();
                fps_convolve(ring, f, g, n)
            },
            0,
            Some(n),
        )
    }

    /// `R[x]`, the polynomials under [`poly_convolve`].
    pub fn poly(ring: R) -> Self {
        Self::new(ring, poly_convolve, 0, None)
    }
}

impl<R: Semiring<Value: Clone>, F: Fn(&R, Vec<R::Value>, Vec<R::Value>) -> Vec<R::Value>>
    MonoidAlgebra<R, F>
{
    fn check(&self, a: &[R::Value]) {
        if let Some(order) = self.order {
            assert!(
                a.len() == order,
                "length mismatch: len={}, order={order}",
                a.len(),
            );
        }
    }
}

impl<R: Semiring<Value: Clone>, F: Fn(&R, Vec<R::Value>, Vec<R::Value>) -> Vec<R::Value>> Semiring
    for MonoidAlgebra<R, F>
{
    type Value = Vec<R::Value>;
    fn zero(&self) -> Vec<R::Value> {
        (0..self.order.unwrap_or(0))
            .map(|_| self.ring.zero())
            .collect()
    }
    fn one(&self) -> Vec<R::Value> {
        let mut e: Vec<R::Value> = (0..self.order.unwrap_or(self.id + 1))
            .map(|_| self.ring.zero())
            .collect();
        e[self.id] = self.ring.one();
        e
    }
    fn add(&self, a: &Vec<R::Value>, b: &Vec<R::Value>) -> Vec<R::Value> {
        self.check(a);
        self.check(b);
        let (long, short) = if a.len() >= b.len() { (a, b) } else { (b, a) };
        let mut c = long.clone();
        for (x, y) in c.iter_mut().zip(short) {
            *x = self.ring.add(x, y);
        }
        c
    }
    fn mul(&self, a: &Vec<R::Value>, b: &Vec<R::Value>) -> Vec<R::Value> {
        self.check(a);
        self.check(b);
        (self.mul)(&self.ring, a.clone(), b.clone())
    }
}
impl<R: Ring<Value: Clone>, F: Fn(&R, Vec<R::Value>, Vec<R::Value>) -> Vec<R::Value>> Ring
    for MonoidAlgebra<R, F>
{
    fn neg(&self, a: &Vec<R::Value>) -> Vec<R::Value> {
        self.check(a);
        a.iter().map(|x| self.ring.neg(x)).collect()
    }
}
