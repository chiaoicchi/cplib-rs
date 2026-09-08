use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::monge::monotone_minima;

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
    let a: Vec<u32> = (0..n).map(|_| parse!(u32)).collect();
    let b: Vec<u32> = (0..m).map(|_| parse!(u32)).collect();
    let idx = monotone_minima(n + m - 1, m, |i: usize, j: usize| {
        if i < j || i - j >= n {
            std::u32::MAX
        } else {
            a[i - j] + b[j]
        }
    });
    for (i, &j) in idx.iter().enumerate() {
        write!(stdout, "{} ", a[i - j] + b[j]).ok();
    }
    writeln!(stdout).ok();
}
