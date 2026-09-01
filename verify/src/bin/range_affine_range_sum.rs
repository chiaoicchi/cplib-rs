use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::algebra::{FnAction, FnMonoid};
use cplib::collections::lazy_segment_tree::LazySegmentTree;

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
    let a: Vec<(u32, u32)> = (0..n).map(|_| (parse!(u32), 1)).collect();
    let mut lazy_segment_tree = LazySegmentTree::from_vec(
        FnMonoid {
            id: (0, 0),
            op: |(a, b): &(u32, u32), (c, d): &(u32, u32)| -> (u32, u32) {
                let x = a + c;
                (if x >= MOD { x - MOD } else { x }, b + d)
            },
        },
        FnMonoid {
            id: (1, 0),
            op: |(a, b): &(u32, u32), (c, d): &(u32, u32)| -> (u32, u32) {
                let f = (*a as u64 * *c as u64 % MOD as u64) as u32;
                let g = (*c as u64 * *b as u64 % MOD as u64) as u32 + d;
                (f, if g >= MOD { g - MOD } else { g })
            },
        },
        FnAction {
            act: |(f, g): &(u32, u32), (a, b): &(u32, u32)| -> (u32, u32) {
                let x = (*f as u64 * *a as u64 % MOD as u64) as u32
                    + (*g as u64 * *b as u64 % MOD as u64) as u32;
                (if x >= MOD { x - MOD } else { x }, *b)
            },
        },
        a,
    );

    for _ in 0..q {
        if parse!(u8) == 0 {
            let l = parse!(usize);
            let r = parse!(usize);
            let c = parse!(u32);
            let d = parse!(u32);
            lazy_segment_tree.range_apply(l..r, &(c, d));
        } else {
            let l = parse!(usize);
            let r = parse!(usize);
            let ans = lazy_segment_tree.fold(l..r).0;
            writeln!(stdout, "{ans}").ok();
        }
    }
}
