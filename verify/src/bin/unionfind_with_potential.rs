use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::algebra::additive::Additive;
use cplib::algebra::canonical::Canonical;
use cplib::collections::potential_dsu::PotentialDsu;
use cplib::num::fp::fp;

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
    let q = parse!(u32);

    let mut potential_dsu = PotentialDsu::new(Additive(Canonical::new()), n);

    for _ in 0..q {
        let t = parse!(u8);
        let u = parse!(usize);
        let v = parse!(usize);
        if t == 0 {
            let p = parse!(u32);
            let b = potential_dsu.unite(v, u, &fp!(p, mod P));
            writeln!(stdout, "{}", if b { 1 } else { 0 }).ok();
        } else {
            let p = potential_dsu.potential(v, u);
            if let Some(p) = p {
                writeln!(stdout, "{p}").ok();
            } else {
                writeln!(stdout, "-1").ok();
            }
        }
    }
}
