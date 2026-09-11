use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::collections::tree::Tree;

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

    let mut e = vec![vec![]; n];
    for _ in 0..n - 1 {
        let a = parse!(usize);
        let b = parse!(usize);
        e[a].push(b);
        e[b].push(a);
    }
    let tree = Tree::from_adjacency(&e, 0);

    for _ in 0..q {
        let s = parse!(usize);
        let t = parse!(usize);
        let i = parse!(usize);
        writeln!(stdout, "{}", tree.jump(s, t, i).unwrap_or(!0) as isize).ok();
    }
}
