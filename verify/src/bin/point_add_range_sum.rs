use std::io::{BufWriter, Read, Write, stdin, stdout};
use std::num::Wrapping;

use cplib::algebra::additive::Additive;
use cplib::algebra::canonical::Canonical;
use cplib::collections::fenwick_tree::FenwickTree;

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

    let a: Vec<Wrapping<u64>> = (0..n).map(|_| Wrapping(parse!(u64))).collect();
    let mut fenwick_tree = FenwickTree::from_vec(Additive(Canonical::new()), a);

    for _ in 0..q {
        let t = parse!(u8);
        if t == 0 {
            let p = parse!(usize);
            let x = parse!(u64);
            fenwick_tree.op_assign(p, &Wrapping(x));
        } else {
            let l = parse!(usize);
            let r = parse!(usize);
            let ans = fenwick_tree.fold(l..r);
            writeln!(stdout, "{ans}").ok();
        }
    }
}
