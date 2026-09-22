/// Rewrites `a` as the next rearrangement in lexicographic order.
///
/// # Definition
/// List the distinct rearrangements of the multiset `a` in lexicographic order. If `a` is not the
/// last one, `a` becomes its successor and `true` is returned; otherwise `a` becomes the first one,
/// that is `a` sorted in non-decreasing order, and `false` is returned.
///
/// # Complexity
/// - Time: O(n), amortized O(1) over a full cycle if the elements are pairwise distinct
/// - Space: O(1)
pub fn next_permutation<T: Ord>(a: &mut [T]) -> bool {
    next_permutation_by(a, T::cmp)
}

/// Rewrites `a` as the next rearrangement in the lexicographic order induced by `compare`.
///
/// # Definition
/// As `next_permutation`, for the sequence of the classes of the elements of `a` under the
/// equivalence `compare(x, y) = Equal`, ordered by `compare`. The elements of a class may be
/// reordered among themselves.
///
/// # Contract
/// `compare` is a total order: `compare(y, x) = compare(x, y).reverse()`, and
/// `compare(x, z) != Greater` whenever `compare(x, y) != Greater` and `compare(y, z) != Greater`.
///
/// # Complexity
/// - Time: O(n) calls of `compare`, amortized O(1) over a full cycle if the elements are pairwise
///   inequivalent
/// - Space: O(1)
pub fn next_permutation_by<T>(
    a: &mut [T],
    mut compare: impl FnMut(&T, &T) -> std::cmp::Ordering,
) -> bool {
    let Some(i) = a
        .windows(2)
        .rposition(|w| compare(&w[0], &w[1]) == std::cmp::Ordering::Less)
    else {
        a.reverse();
        return false;
    };
    let j = a
        .iter()
        .rposition(|x| compare(x, &a[i]) == std::cmp::Ordering::Greater)
        .unwrap();
    a.swap(i, j);
    a[i + 1..].reverse();
    true
}

/// Rewrites `a` as the next rearrangement in the lexicographic order of the keys `f`.
///
/// # Definition
/// As `next_permutation_by` with `compare(x, y) = f(x).cmp(&f(y))`.
///
/// # Complexity
/// - Time: O(n) calls of `f`, amortized O(1) over a full cycle if the keys are pairwise distinct
/// - Space: O(1)
pub fn next_permutation_by_key<T, K: Ord>(a: &mut [T], mut f: impl FnMut(&T) -> K) -> bool {
    next_permutation_by(a, |x, y| f(x).cmp(&f(y)))
}

/// Rewrites `a` as the previous rearrangement in lexicographic order.
///
/// # Definition
/// List the distinct rearrangements of the multiset `a` in lexicographic order. If `a` is not the
/// first one, `a` becomes its predecessor and `true` is returned; otherwise `a` becomes the last
/// one, that is `a` sorted in non-increasing order, and `false` is returned.
///
/// # Complexity
/// - Time: O(n), amortized O(1) over a full cycle if the elements are pairwise distinct
/// - Space: O(1)
pub fn prev_permutation<T: Ord>(a: &mut [T]) -> bool {
    prev_permutation_by(a, T::cmp)
}

/// Rewrites `a` as the previous rearrangement in the lexicographic order induced by `compare`.
///
/// # Definition
/// As `prev_permutation`, for the sequence of the classes of the elements of `a` under the
/// equivalence `compare(x, y) = Equal`, ordered by `compare`. The elements of a class may be
/// reordered among themselves. Equivalently, it is `next_permutation_by` for the opposite order.
///
/// # Contract
/// `compare` is a total preorder, as in `next_permutation_by`.
///
/// # Complexity
/// - Time: O(n) calls of `compare`, amortized O(1) over a full cycle if the elements are pairwise
///   inequivalent
/// - Space: O(1)
pub fn prev_permutation_by<T>(
    a: &mut [T],
    mut compare: impl FnMut(&T, &T) -> std::cmp::Ordering,
) -> bool {
    next_permutation_by(a, |x, y| compare(y, x))
}

/// Rewrites `a` as the previous rearrangement in the lexicographic order of the keys `f`.
///
/// # Definition
/// As `prev_permutation_by` with `compare(x, y) = f(x).cmp(&f(y))`.
///
/// # Complexity
/// - Time: O(n) calls of `f`, amortized O(1) over a full cycle if the keys are pairwise distinct
/// - Space: O(1)
pub fn prev_permutation_by_key<T, K: Ord>(a: &mut [T], mut f: impl FnMut(&T) -> K) -> bool {
    prev_permutation_by(a, |x, y| f(x).cmp(&f(y)))
}
