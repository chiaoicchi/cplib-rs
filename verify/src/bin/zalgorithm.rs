use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::sequence::z_algorithm::z_algorithm;

fn main() {
    let mut input = Vec::new();
    stdin().lock().read_to_end(&mut input).unwrap();
    let mut iter = input.split(|&b| b <= b' ').filter(|s| !s.is_empty());
    let mut stdout = BufWriter::new(stdout().lock());

    let s: &[u8] = iter.next().unwrap();
    let z = z_algorithm(s);
    for z in z {
        write!(stdout, "{z} ").ok();
    }
    writeln!(stdout).ok();
}
