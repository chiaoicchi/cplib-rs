use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::string::suffix_array::SuffixArray;

fn main() {
    let mut input = Vec::new();
    stdin().lock().read_to_end(&mut input).unwrap();
    let mut iter = input.split(|&b| b <= b' ').filter(|s| !s.is_empty());
    let mut stdout = BufWriter::new(stdout().lock());

    let s: &[u8] = iter.next().unwrap();
    let sa = SuffixArray::from_slice(s);
    for a in sa.sa() {
        write!(stdout, "{a} ").ok();
    }
    writeln!(stdout).ok();
}
