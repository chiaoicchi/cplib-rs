use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::sequence::manacher::manacher;

fn main() {
    let mut input = Vec::new();
    stdin().lock().read_to_end(&mut input).unwrap();
    let mut iter = input.split(|&b| b <= b' ').filter(|s| !s.is_empty());
    let mut stdout = BufWriter::new(stdout().lock());

    let s: &[u8] = iter.next().unwrap();
    let r = manacher(
        &s.iter()
            .flat_map(|&c| [b'A', c])
            .chain(b"A".iter().copied())
            .collect::<Vec<_>>(),
    );
    for r in &r[1..r.len() - 1] {
        write!(stdout, "{} ", r - 1).ok();
    }
    writeln!(stdout).ok();
}
