use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::algebra::affine::Affine;
use cplib::collections::segment_tree::SegmentTree;
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

    let a: Vec<(Fp<P>, Fp<P>)> = (0..n)
        .map(|_| (Fp::new(parse!(u32)), Fp::new(parse!(u32))))
        .collect();
    let mut segment_tree = SegmentTree::from_vec(Affine::new(), a);

    for _ in 0..q {
        let t = parse!(u8);
        if t == 0 {
            let p = parse!(usize);
            let c = parse!(u32);
            let d = parse!(u32);
            segment_tree.set(p, (Fp::new(c), Fp::new(d)));
        } else {
            let l = parse!(usize);
            let r = parse!(usize);
            let x = Fp::new(parse!(u32));
            let (c, d) = segment_tree.fold(l..r);
            let ans = c * x + d;
            writeln!(stdout, "{ans}").ok();
        }
    }
}
