use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::algebra::Monoid;
use cplib::algebra::affine::Affine;
use cplib::algebra::canonical::Canonical;
use cplib::collections::segment_tree::SegmentTree;
use cplib::collections::tree::{Segment, Tree};
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

    let ab: Vec<(Fp<P>, Fp<P>)> = (0..n)
        .map(|_| (fp!(parse!(u32)), fp!(parse!(u32))))
        .collect();
    let mut e = vec![vec![]; n];
    for _ in 0..n - 1 {
        let u = parse!(usize);
        let v = parse!(usize);
        e[u].push(v);
        e[v].push(u);
    }
    let tree = Tree::from_adjacency(&e, 0);
    let monoid = Affine(Canonical::new());
    let mut segment_tree =
        SegmentTree::from_vec(monoid, (0..n).map(|i| ab[tree.vertex(i)]).collect());
    let mut inv_segment_tree =
        SegmentTree::from_vec(monoid, (0..n).rev().map(|i| ab[tree.vertex(i)]).collect());

    for _ in 0..q {
        let t = parse!(u8);
        match t {
            0 => {
                let p = parse!(usize);
                let c = parse!(u32);
                let d = parse!(u32);
                segment_tree.set(tree.index(p), (fp!(c), fp!(d)));
                inv_segment_tree.set(n - 1 - tree.index(p), (fp!(c), fp!(d)));
            }
            1 => {
                let u = parse!(usize);
                let v = parse!(usize);
                let x = parse!(u32);
                let mut acc = monoid.id();
                for seg in tree.path(u, v) {
                    let x = match seg {
                        Segment::Up(r) => inv_segment_tree.fold(n - r.end..n - r.start),
                        Segment::Down(r) => segment_tree.fold(r),
                    };
                    acc = monoid.op(&acc, &x);
                }
                writeln!(stdout, "{}", acc.0 * fp!(x) + acc.1).ok();
            }
            _ => unreachable!(),
        }
    }
}
