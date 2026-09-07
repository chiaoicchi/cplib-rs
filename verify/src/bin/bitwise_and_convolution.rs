use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::algebra::and::And;
use cplib::algebra::canonical::Canonical;
use cplib::convolution::convolve_by_transform;
use cplib::num::fp::Fp;

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
    let a: Vec<Fp<P>> = (0..1 << n).map(|_| Fp::new(parse!(u32))).collect();
    let b: Vec<Fp<P>> = (0..1 << n).map(|_| Fp::new(parse!(u32))).collect();

    let ans = convolve_by_transform(&And::new(), &Canonical::new(), a, b);
    for ans in ans {
        write!(stdout, "{ans} ").ok();
    }
    writeln!(stdout).ok();
}
