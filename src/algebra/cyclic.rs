use crate::algebra::Monoid;

/// The additive group `Z/nZ`.
///
/// # Definition
/// `(Z/nZ, 0, +)` is the cyclic group of order `n`. Identifying `[0, n)` with the residues, the
/// operation is `(a + b) % n`.
///
/// # Contract
/// `n >= 1`
#[derive(Clone, Copy)]
pub struct Cyclic {
    n: usize,
}
impl Cyclic {
    pub const fn new(n: usize) -> Self {
        Self { n }
    }
    #[allow(clippy::len_without_is_empty)]
    pub const fn len(&self) -> usize {
        self.n
    }
}

impl Monoid for Cyclic {
    type Value = usize;
    fn id(&self) -> usize {
        0
    }
    fn op(&self, a: &usize, b: &usize) -> usize {
        (a + b) % self.n
    }
}
