use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::graph::min_cost_flow::MinCostFlow;

fn main() {
    let mut input = Vec::new();
    stdin().lock().read_to_end(&mut input).unwrap();
    let mut iter = input.split(|&b| b <= b' ').filter(|s| !s.is_empty());
    let mut stdout = BufWriter::new(stdout().lock());

    macro_rules! parse {
        ($t:ty) => {{
            let s = iter.next().unwrap();
            let (neg, digits) = match s[0] {
                b'-' => (true, &s[1..]),
                _ => (false, s),
            };
            let mut x: $t = 0;
            for &b in digits {
                x = x * 10 + (b - b'0') as $t;
            }
            if neg { 0 - x } else { x }
        }};
    }

    let n = parse!(usize);
    let mut network = MinCostFlow::new(2 * n + 2);
    for i in 0..n {
        for j in 0..n {
            let a = parse!(i64);
            network.add_edge(i, n + j, 1, a);
        }
    }
    for i in 0..n {
        network.add_edge(2 * n, i, 1u32, 0i64);
        network.add_edge(n + i, 2 * n + 1, 1, 0);
    }
    let cost = network.flow_exact(2 * n, 2 * n + 1, n as u32).unwrap();
    writeln!(stdout, "{cost}").ok();
    for e in 0..n * n {
        let (_, j, _, f, _) = network.get_edge(e);
        if f == 1 {
            write!(stdout, "{} ", j - n).ok();
        }
    }
}
