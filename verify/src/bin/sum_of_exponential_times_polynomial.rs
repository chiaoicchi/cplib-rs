use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::algebra::canonical::{Canonical, Multiplicative};
use cplib::arithmetic::multiplicative::multiplicative_table;
use cplib::num::fp::{Fp, fp};
use cplib::poly::geometric::iota_geometric_sum;

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

    let r = parse!(u32);
    let d = parse!(usize);
    let n = parse!(u64);
    let mut y = multiplicative_table(&Multiplicative(Canonical::<Fp<P>>::new()), d, |p, e| {
        fp!(p.pow(e) as u32).pow(d as u64)
    });
    y[0] = fp!((d == 0) as u32);
    let ans = iota_geometric_sum(&Canonical::new(), &y, &fp!(r), n);
    writeln!(stdout, "{ans}").ok();
}
