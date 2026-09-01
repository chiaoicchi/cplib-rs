use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::algebra::FnMonoid;
use cplib::collections::segment_tree::SegmentTree;

const MOD: u32 = 998_244_353;

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

    let a: Vec<(u32, u32)> = (0..n).map(|_| (parse!(u32), parse!(u32))).collect();
    let mut segment_tree = SegmentTree::from_vec(
        FnMonoid {
            id: (1, 0),
            op: |&(a, b): &(u32, u32), &(c, d): &(u32, u32)| -> (u32, u32) {
                let f = (a as u64 * c as u64 % MOD as u64) as u32;
                let g = (c as u64 * b as u64 % MOD as u64) as u32 + d;
                (f, if g >= MOD { g - MOD } else { g })
            },
        },
        a,
    );

    for _ in 0..q {
        let t = parse!(u8);
        if t == 0 {
            let p = parse!(usize);
            let c = parse!(u32);
            let d = parse!(u32);
            segment_tree.set(p, (c, d));
        } else {
            let l = parse!(usize);
            let r = parse!(usize);
            let x = parse!(u64);
            let (c, d) = segment_tree.fold(l..r);
            let ans = (c as u64 * x % MOD as u64) as u32 + d;
            writeln!(stdout, "{}", if ans >= MOD { ans - MOD } else { ans }).ok();
        }
    }
}
