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

    let n = parse!(usize);
    let q = parse!(u32);
    let a: Vec<u64> = (0..n).map(|_| parse!(u64)).collect();
    let mut binary_trie = BinaryTrie::from_vec(30, a);

    for _ in 0..q {
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
                let ans = binary_trie.kth(x - 1).unwrap_or(!0);
                writeln!(stdout, "{}", ans as i64).ok();
            }
            3 => {
                let ans = binary_trie.rank(x + 1);
                writeln!(stdout, "{}", ans).ok();
            }
            4 => {
                let ans = binary_trie.floor(x).unwrap_or(!0);
                writeln!(stdout, "{}", ans as i64).ok();
            }
            5 => {
                let ans = binary_trie.ceil(x).unwrap_or(!0);
                writeln!(stdout, "{}", ans as i64).ok();
            }
            _ => unreachable!(),
        }
    }
}
