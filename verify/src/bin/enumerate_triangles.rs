use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::graph::short_cycle::for_each_triangle;
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
    let m = parse!(usize);
    let x: Vec<Fp<P>> = (0..n).map(|_| fp!(parse!(u32), mod P)).collect();
    let e: Vec<(usize, usize)> = (0..m).map(|_| (parse!(usize), parse!(usize))).collect();
    let mut ans = fp!(0);
    for_each_triangle(n, &e, |a: usize, b: usize, c: usize| {
        ans += x[a] * x[b] * x[c]
    });
    writeln!(stdout, "{ans}").ok();
}
