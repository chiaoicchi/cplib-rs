use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::collections::suffix_array::SuffixArray;

fn main() {
    let mut input = Vec::new();
    stdin().lock().read_to_end(&mut input).unwrap();
    let mut iter = input.split(|&b| b <= b' ').filter(|s| !s.is_empty());
    let mut stdout = BufWriter::new(stdout().lock());

    let s: &[u8] = iter.next().unwrap();
    let sa = SuffixArray::from_slice(s);
    let ans = s.len() * (s.len() + 1) / 2 - sa.lcp().iter().sum::<usize>();
    writeln!(stdout, "{ans}").ok();
}
