use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::collections::suffix_array::SuffixArray;

fn main() {
    let mut input = Vec::new();
    stdin().lock().read_to_end(&mut input).unwrap();
    let mut iter = input.split(|&b| b <= b' ').filter(|s| !s.is_empty());
    let mut stdout = BufWriter::new(stdout().lock());

    let s: &[u8] = iter.next().unwrap();
    let n = s.len();
    let t: &[u8] = iter.next().unwrap();
    let sa = SuffixArray::from_slice(&s.iter().chain(b"A").chain(t).copied().collect::<Vec<_>>());
    let mut ans_val = 0;
    let mut ans_pos = (0, 0, 0, 0);
    for (v, &lcp) in sa.sa().windows(2).zip(sa.lcp().iter()) {
        let i = v[0];
        let j = v[1];
        let (a, b) = if i < n { (i, j) } else { (j, i) };
        if a < n && n < b && ans_val < lcp {
            ans_val = lcp;
            ans_pos = (a, a + lcp, b - n - 1, b - n - 1 + lcp);
        }
    }
    writeln!(
        stdout,
        "{} {} {} {}",
        ans_pos.0, ans_pos.1, ans_pos.2, ans_pos.3
    )
    .ok();
}
