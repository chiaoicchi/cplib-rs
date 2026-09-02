use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::algebra::closures::{FnAction, FnMonoid};
use cplib::collections::dual_segment_tree::DualSegmentTree;

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

    let a: Vec<u32> = (0..n).map(|_| parse!(u32)).collect();
    let mut dual_segment_tree = DualSegmentTree::from_vec(
        FnMonoid {
            id: (1, 0),
            op: |&(a, b): &(u32, u32), &(c, d): &(u32, u32)| -> (u32, u32) {
                let f = (a as u64 * c as u64 % MOD as u64) as u32;
                let g = (c as u64 * b as u64 % MOD as u64) as u32 + d;
                (f, if g >= MOD { g - MOD } else { g })
            },
        },
        FnAction {
            act: |(a, b): &(u32, u32), x: &u32| -> u32 {
                let x = (*a as u64 * *x as u64 % MOD as u64) as u32 + b;
                if x >= MOD { x - MOD } else { x }
            },
        },
        a,
    );

    for _ in 0..q {
        let t = parse!(u8);
        if t == 0 {
            let l = parse!(usize);
            let r = parse!(usize);
            let b = parse!(u32);
            let c = parse!(u32);
            dual_segment_tree.range_apply(l..r, &(b, c));
        } else {
            let i = parse!(usize);
            let ans = dual_segment_tree.get(i).unwrap();
            writeln!(stdout, "{ans}").ok();
        }
    }
}
