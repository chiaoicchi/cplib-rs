use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::collections::binary_trie::BinaryTrie;

fn main() {
    let mut input = Vec::new();
    stdin().lock().read_to_end(&mut input).unwrap();
    let mut iter = input.split(|&b| b <= b' ').filter(|s| !s.is_empty());
    let mut stdout = BufWriter::new(stdout().lock());

    macro_rules! parse {
        ($t:ty) => {{
            let s = iter.next().unwrap();
            let mut x: $t = 0;
            for &b in s {
                x = x * 10 + (b - b'0') as $t;
            }
            x
        }};
    }

    let mut binary_trie = BinaryTrie::new(30);

    for _ in 0..parse!(u32) {
        let c = parse!(u8);
        let x = parse!(u64);
        match c {
            0 => {
                if !binary_trie.contains(x) {
                    binary_trie.insert(x);
                }
            }
            1 => {
                binary_trie.remove(x);
            }
            2 => {
                let ans = binary_trie.min_xor(x).unwrap();
                writeln!(stdout, "{}", ans).ok();
            }
            _ => unreachable!(),
        }
    }
}
