use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::algebra::affine::Affine;
use cplib::algebra::closures::{FnAction, FnMonoid};
use cplib::collections::lazy_segment_tree::LazySegmentTree;
use cplib::num::fp::Fp;

const P: u32 = 998_244_353;

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
    let a: Vec<(Fp<P>, Fp<P>)> = (0..n).map(|_| (Fp::new(parse!(u32)), Fp::new(1))).collect();
    let mut lazy_segment_tree = LazySegmentTree::from_vec(
        FnMonoid {
            id: (Fp::new(0), Fp::new(0)),
            op: |(a, b): &(Fp<P>, Fp<P>), (c, d): &(Fp<P>, Fp<P>)| -> (Fp<P>, Fp<P>) {
                (a + c, b + d)
            },
        },
        Affine::new(),
        FnAction {
            act: |(f, g): &(Fp<P>, Fp<P>), (a, b): &(Fp<P>, Fp<P>)| -> (Fp<P>, Fp<P>) {
                (f * a + g * b, *b)
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
            lazy_segment_tree.range_apply(l..r, &(Fp::new(c), Fp::new(d)));
        } else {
            let l = parse!(usize);
            let r = parse!(usize);
            let ans = lazy_segment_tree.fold(l..r).0;
            writeln!(stdout, "{ans}").ok();
        }
    }
}
