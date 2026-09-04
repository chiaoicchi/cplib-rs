use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::algebra::canonical::Canonical;
use cplib::linear::dyn_matrix::DynMatrix;
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
    let m = parse!(usize);
    let k = parse!(usize);
    let a: Vec<Vec<Fp<P>>> = (0..n)
        .map(|_| (0..m).map(|_| Fp::new(parse!(u32))).collect())
        .collect();
    let b: Vec<Vec<Fp<P>>> = (0..m)
        .map(|_| (0..k).map(|_| Fp::new(parse!(u32))).collect())
        .collect();

    let mat_a = Matrix::from_vec(Canonical::new(), a);
    let mat_b = Matrix::from_vec(Canonical::new(), b);

    let ans = mat_a * mat_b;
    for vi in ans.iter() {
        for vij in vi.iter() {
            write!(stdout, "{vij} ").ok();
        }
        writeln!(stdout).ok();
    }
}
