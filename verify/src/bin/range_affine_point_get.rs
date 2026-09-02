use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::algebra::affine::Affine;
use cplib::algebra::canonical::Canonical;
use cplib::collections::dual_segment_tree::DualSegmentTree;
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

    let a: Vec<Fp<P>> = (0..n).map(|_| Fp::new(parse!(u32))).collect();
    let mut dual_segment_tree =
        DualSegmentTree::from_vec(Affine(Canonical::new()), Affine(Canonical::new()), a);

    for _ in 0..q {
        let t = parse!(u8);
        if t == 0 {
            let l = parse!(usize);
            let r = parse!(usize);
            let b = parse!(u32);
            let c = parse!(u32);
            dual_segment_tree.range_apply(l..r, &(Fp::new(b), Fp::new(c)));
        } else {
            let i = parse!(usize);
            let ans = dual_segment_tree.get(i).unwrap();
            writeln!(stdout, "{ans}").ok();
        }
    }
}
