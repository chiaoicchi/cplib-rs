use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::graph::short_cycle::count_c4;

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
    let e: Vec<(usize, usize)> = (0..m).map(|_| (parse!(usize), parse!(usize))).collect();
    let ans = count_c4(n, &e);
    for ans in ans {
        write!(stdout, "{ans} ").ok();
    }
    writeln!(stdout).ok();
}
