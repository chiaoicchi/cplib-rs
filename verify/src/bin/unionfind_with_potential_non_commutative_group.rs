use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::algebra::closures::FnGroup;
use cplib::collections::potential_dsu::PotentialDsu;
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
    let q = parse!(u32);

    let mut potential_dsu = PotentialDsu::new(
        FnGroup {
            id: [fp!(1), fp!(0), fp!(0), fp!(1)],
            op: |a: &[Fp<P>; 4], b: &[Fp<P>; 4]| -> [Fp<P>; 4] {
                [
                    a[0] * b[0] + a[1] * b[2],
                    a[0] * b[1] + a[1] * b[3],
                    a[2] * b[0] + a[3] * b[2],
                    a[2] * b[1] + a[3] * b[3],
                ]
            },
            inv: |a: &[Fp<P>; 4]| -> [Fp<P>; 4] { [a[3], -a[1], -a[2], a[0]] },
        },
        n,
    );

    for _ in 0..q {
        let t = parse!(u8);
        let u = parse!(usize);
        let v = parse!(usize);
        if t == 0 {
            let x = [
                fp!(parse!(u32)),
                fp!(parse!(u32)),
                fp!(parse!(u32)),
                fp!(parse!(u32)),
            ];
            let b = potential_dsu.unite(v, u, &x);
            writeln!(stdout, "{}", if b { 1 } else { 0 }).ok();
        } else {
            let p = potential_dsu.potential(v, u);
            if let Some(p) = p {
                writeln!(stdout, "{} {} {} {}", p[0], p[1], p[2], p[3]).ok();
            } else {
                writeln!(stdout, "-1").ok();
            }
        }
    }
}
