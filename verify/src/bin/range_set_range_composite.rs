use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::algebra::canonical::{Affine, Canonical};
use cplib::algebra::closures::{FnAction, FnMonoid};
use cplib::algebra::power::pow;
use cplib::algebra::{Monoid, Semigroup};
use cplib::collections::lazy_segment_tree::LazySegmentTree;
use cplib::num::fp::{Fp, fp};

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
    let affine = Affine(Canonical::<Fp<P>>::new());
    let mut lazy_segment_tree = LazySegmentTree::from_vec(
        FnMonoid {
            id: (affine.id(), 0),
            op: |x: &((Fp<P>, Fp<P>), usize),
                 y: &((Fp<P>, Fp<P>), usize)|
             -> ((Fp<P>, Fp<P>), usize) { (affine.op(&x.0, &y.0), x.1 + y.1) },
        },
        FnMonoid {
            id: None,
            op: |a: &Option<(Fp<P>, Fp<P>)>,
                 b: &Option<(Fp<P>, Fp<P>)>|
             -> Option<(Fp<P>, Fp<P>)> { if b.is_some() { *b } else { *a } },
        },
        FnAction {
            act: |f: &Option<(Fp<P>, Fp<P>)>,
                  x: &((Fp<P>, Fp<P>), usize)|
             -> ((Fp<P>, Fp<P>), usize) {
                match f {
                    Some(g) => (pow(&affine, g, x.1 as u64), x.1),
                    None => *x,
                }
            },
        },
        a.into_iter().map(|f| (f, 1)).collect(),
    );

    for _ in 0..q {
        if parse!(u8) == 0 {
            let l = parse!(usize);
            let r = parse!(usize);
            let c = parse!(u32);
            let d = parse!(u32);
            lazy_segment_tree.range_apply(l..r, &Some((Fp::new(c), Fp::new(d))));
        } else {
            let l = parse!(usize);
            let r = parse!(usize);
            let x = parse!(u32);
            let ((c, d), _) = lazy_segment_tree.fold(l..r);
            writeln!(stdout, "{}", c * fp!(x) + d).ok();
        }
    }
}
