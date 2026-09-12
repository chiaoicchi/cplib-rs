use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::algebra::canonical::Canonical;
use cplib::num::fp::Fp;
use cplib::set_power_series::SetPowerSeries;

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
    let b: Vec<Fp<P>> = (0..1 << n).map(|_| Fp::new(parse!(u32))).collect();

    let ans = SetPowerSeries::new(Canonical::new(), n).exp(&b);
    for ans in ans {
        write!(stdout, "{ans} ").ok();
    }
    writeln!(stdout).ok();
}
