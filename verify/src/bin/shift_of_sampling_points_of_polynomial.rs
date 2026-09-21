use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::algebra::canonical::Canonical;
use cplib::num::fp::{Fp, fp};
use cplib::poly::interpolation::shift_of_sampling_points;

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
    let c = Fp::<P>::new(parse!(u32));
    let f: Vec<Fp<P>> = (0..n).map(|_| fp!(parse!(u32))).collect();
    let ans = shift_of_sampling_points(&Canonical::new(), &f, &c, m);
    for ans in ans {
        write!(stdout, "{ans} ").ok();
    }
    writeln!(stdout).ok();
}
