use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::graph::spanning_tree::minimum_spanning_tree;

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
    let e: Vec<(usize, usize, u64)> = (0..m)
        .map(|_| (parse!(usize), parse!(usize), parse!(u64)))
        .collect();

    let mst = minimum_spanning_tree(n, &e);
    let ans: u64 = (0..m).map(|i| if mst[i] { e[i].2 } else { 0 }).sum();
    writeln!(stdout, "{ans}").ok();
    for (i, b) in mst.iter().enumerate() {
        if *b {
            write!(stdout, "{i} ").ok();
        }
    }
    writeln!(stdout).ok();
}
