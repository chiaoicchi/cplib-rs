/// The half-open interval `[l, r)` given by `range` within `[0, len)`.
///
/// # Definition
/// `l` is `0` for an unbounded start, and `r` is `len` for an unbounded end.
///
/// # Complexity
/// - Time: O(1)
/// - Space: O(1)
///
/// # Panics
/// Panics if `l > r` or `r > len`.
pub(crate) fn to_half_open(len: usize, range: impl std::ops::RangeBounds<usize>) -> (usize, usize) {
    use std::ops::Bound;
    let l = match range.start_bound() {
        Bound::Unbounded => 0,
        Bound::Included(&x) => x,
        Bound::Excluded(&x) => x.saturating_add(1),
    };
    let r = match range.end_bound() {
        Bound::Unbounded => len,
        Bound::Included(&x) => x.saturating_add(1),
        Bound::Excluded(&x) => x,
    };
    assert!(
        l <= r,
        "left bound must be less than or equal to right bound: l={l}, r={r}"
    );
    assert!(r <= len, "range out of bounds: range=[{l}, {r})");
    (l, r)
}
