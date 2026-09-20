/// The shape of a grid with `h` rows and `w` columns.
///
/// # Definition
/// Represents the set of cells `[0, h) x [0, w)`. A cell is a pair `(i, j)` of a row `i` and a
/// column `j`. The contents of the grid are not stored.
///
/// An offset `(di, dj)` is a pair of `usize` interpreted modulo `2^usize::BITS`, so that `!0`
/// stands for `-1`, `!1` for `-2`, and so on.
#[derive(Clone, Copy)]
pub struct GridShape {
    h: usize,
    w: usize,
}

impl GridShape {
    const DIR4: [(usize, usize); 4] = [(!0, 0), (0, !0), (0, 1), (1, 0)];
    const DIR8: [(usize, usize); 8] = [
        (!0, !0),
        (!0, 0),
        (!0, 1),
        (0, !0),
        (0, 1),
        (1, !0),
        (1, 0),
        (1, 1),
    ];

    /// Constructs the shape of a grid with `h` rows and `w` columns.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `h * w` overflows `usize`.
    /// Panics if `h` or `w` is greater than `2^(usize::BITS - 1)`.
    pub fn new(h: usize, w: usize) -> Self {
        assert!(h.checked_mul(w).is_some(), "h * w overflows: h={h}, w={w}");
        let half = 1 << (usize::BITS - 1);
        assert!(
            h <= half && w <= half,
            "h and w must be at most 2^(usize::BITS - 1): h={h}, w={w}"
        );
        Self { h, w }
    }

    /// Returns the number of rows.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn h(&self) -> usize {
        self.h
    }

    /// Returns the number of columns.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    pub fn w(&self) -> usize {
        self.w
    }

    /// Returns the row-major index `i * w + j` of the cell `(i, j)`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i >= h` or `j >= w`.
    pub fn index(&self, i: usize, j: usize) -> usize {
        assert!(
            i < self.h && j < self.w,
            "cell out of bounds: (i,j)=({i},{j}), (h,w)=({},{})",
            self.h,
            self.w
        );
        i * self.w + j
    }

    /// Returns the cell `(u / w, u % w)` whose row-major index is `u`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `u >= h * w`.
    pub fn coord(&self, u: usize) -> (usize, usize) {
        assert!(
            u < self.h * self.w,
            "index out of bounds: u={u}, h={}, w={}",
            self.h,
            self.w
        );
        (u / self.w, u % self.w)
    }

    /// Returns the cell `(i + di, j + dj)`, or `None` if it is outside the grid.
    ///
    /// # Contract
    /// The offsets `di` and `dj`, as signed integers, have absolute value at most
    /// `2^(usize::BITS - 1)`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i >= h` or `j >= w`.
    pub fn shift(
        &self,
        (i, j): (usize, usize),
        (di, dj): (usize, usize),
    ) -> Option<(usize, usize)> {
        assert!(
            i < self.h && j < self.w,
            "cell out of bounds: (i,j)=({i},{j}), (h,w)=({},{})",
            self.h,
            self.w
        );
        let (ni, nj) = (i.wrapping_add(di), j.wrapping_add(dj));
        (ni < self.h && nj < self.w).then_some((ni, nj))
    }

    /// Returns the 4-neighbors of the cell `(i, j)`.
    ///
    /// # Definition
    /// Yields the cells `(i', j')` of the grid with `|i - i'| + |j - j'| = 1`, in the
    /// lexicographic order: `(i - 1, j)`, `(i, j - 1)`, `(i, j + 1)`, `(i + 1, j)`.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i >= h` or `j >= w`.
    pub fn neighbors4(&self, i: usize, j: usize) -> impl Iterator<Item = (usize, usize)> {
        assert!(
            i < self.h && j < self.w,
            "cell out of bounds: (i,j)=({i},{j}), (h,w)=({},{})",
            self.h,
            self.w
        );
        Self::DIR4
            .into_iter()
            .filter_map(move |d| self.shift((i, j), d))
    }

    /// Returns the 8-neighbors of the cell `(i, j)`.
    ///
    /// # Definition
    /// Yields the cells `(i', j')` of the grid with `max(|i - i'|, |j - j'|) = 1`, in the
    /// lexicographic order.
    ///
    /// # Complexity
    /// - Time: O(1)
    /// - Space: O(1)
    ///
    /// # Panics
    /// Panics if `i >= h` or `j >= w`.
    pub fn neighbors8(&self, i: usize, j: usize) -> impl Iterator<Item = (usize, usize)> {
        assert!(
            i < self.h && j < self.w,
            "cell out of bounds: (i,j)=({i},{j}), (h,w)=({},{})",
            self.h,
            self.w
        );
        Self::DIR8
            .into_iter()
            .filter_map(move |d| self.shift((i, j), d))
    }
}
