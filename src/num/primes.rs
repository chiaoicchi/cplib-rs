/// The primes below `n` by the sieve of Eratosthenes.
///
/// # Complexity
/// - Time: O(n log log n)
/// - Space: O(n)
pub fn primes(n: usize) -> Vec<usize> {
    let mut is_prime = vec![true; n];
    let mut ps = Vec::new();
    for p in 2..n {
        if !is_prime[p] {
            continue;
        }
        ps.push(p);
        let mut m = p * p;
        while m < n {
            is_prime[m] = false;
            m += p;
        }
    }
    ps
}
