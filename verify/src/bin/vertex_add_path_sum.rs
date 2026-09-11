use std::io::{BufWriter, Read, Write, stdin, stdout};
use std::num::Wrapping;

use cplib::algebra::additive::Additive;
use cplib::algebra::canonical::Canonical;
use cplib::collections::fenwick_tree::FenwickTree;
use cplib::collections::tree::{Segment, Tree};

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
    let mut e = vec![vec![]; n];
    for _ in 0..n - 1 {
        let u = parse!(usize);
        let v = parse!(usize);
        e[u].push(v);
        e[v].push(u);
    }
    let tree = Tree::from_adjacency(&e, 0);
    let mut fenwick_tree = FenwickTree::from_vec(
        Additive(Canonical::new()),
        (0..n).map(|i| a[tree.vertex(i)]).collect(),
    );

    for _ in 0..q {
        let t = parse!(u8);
        match t {
            0 => {
                let p = parse!(usize);
                let x = parse!(u64);
                fenwick_tree.op_assign(tree.index(p), &Wrapping(x));
            }
            1 => {
                let u = parse!(usize);
                let v = parse!(usize);
                let mut ans = Wrapping(0);
                for seg in tree.path(u, v) {
                    let (Segment::Up(r) | Segment::Down(r)) = seg;
                    ans += fenwick_tree.fold(r);
                }
                writeln!(stdout, "{ans}").ok();
            }
            _ => unreachable!(),
        }
    }
}
