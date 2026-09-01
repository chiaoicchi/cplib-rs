/// Converts `range` to a pair `(l, r)` representing the half-open interval `[l, r)`,
/// where an unbounded end defaults to `0` or `max`.
pub(crate) fn to_half_open(max: usize, range: impl std::ops::RangeBounds<usize>) -> (usize, usize) {
    use std::ops::Bound;
    let l = match range.start_bound() {
        Bound::Unbounded => 0,
        Bound::Included(&x) => x,
        Bound::Excluded(&x) => x + 1,
    };
    let r = match range.end_bound() {
        Bound::Unbounded => max,
        Bound::Included(&x) => x + 1,
        Bound::Excluded(&x) => x,
    };
    (l, r)
}
