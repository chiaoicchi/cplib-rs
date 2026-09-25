use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::graph::shortest_path::ShortestPathTree;

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
    let m = parse!(usize);
    let s = parse!(usize);
    let t = parse!(usize);
    let mut e = vec![vec![]; n];
    for _ in 0..m {
        let a = parse!(usize);
        let b = parse!(usize);
        let c = parse!(u64);
        e[a].push((b, c));
    }
    let spt = ShortestPathTree::dijkstra(n, &[s], |v: usize| &e[v]);
    let Some(d) = spt.dist(t) else {
        writeln!(stdout, "-1").ok();
        return;
    };
    let path = spt.path(t).unwrap();
    writeln!(stdout, "{} {}", d, path.len() - 1).ok();
    for v in path.windows(2) {
        writeln!(stdout, "{} {}", v[0], v[1]).ok();
    }
}
