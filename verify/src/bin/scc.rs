use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::graph::scc::Scc;

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
    let mut e = vec![vec![]; n];
    for _ in 0..m {
        let a = parse!(usize);
        let b = parse!(usize);
        e[a].push(b);
    }
    let scc = Scc::from_adjacency(&e);
    writeln!(stdout, "{}", scc.len()).ok();
    for i in 0..scc.len() {
        write!(stdout, "{} ", scc.group(i).len()).ok();
        for v in scc.group(i) {
            write!(stdout, "{} ", v).ok();
        }
    }
}
