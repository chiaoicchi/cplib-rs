use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::graph::max_flow::MaxFlow;

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

    let l = parse!(usize);
    let r = parse!(usize);
    let m = parse!(usize);

    let mut max_flow = MaxFlow::new(l + r + 2);

    for _ in 0..m {
        let a = parse!(usize);
        let b = parse!(usize);
        max_flow.add_edge(a, b + l, 1);
    }
    for i in 0..l {
        max_flow.add_edge(l + r, i, 1);
    }
    for i in 0..r {
        max_flow.add_edge(i + l, l + r + 1, 1);
    }
    let ans = max_flow.flow(l + r, l + r + 1);
    writeln!(stdout, "{ans}").ok();
    for i in 0..m {
        let (from, to, _, flow) = max_flow.get_edge(i);
        if flow == 1 {
            writeln!(stdout, "{from} {}", to - l).ok();
        }
    }
}
