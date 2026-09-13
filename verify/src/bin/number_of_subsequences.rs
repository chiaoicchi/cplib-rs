use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::collections::compression::Compression;
use cplib::num::fp::Fp;
use cplib::sequence::num_subsequences::num_subsequences;

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

    let a: Vec<u32> = (0..n).map(|_| parse!(u32)).collect();
    let comp_a = Compression::from_vec(a.clone());
    let ans: Fp<P> = num_subsequences(
        &(0..n)
            .map(|i| comp_a.compress(&a[i]))
            .map(Option::unwrap)
            .collect::<Vec<_>>(),
        n,
    );
    writeln!(stdout, "{}", ans - Fp::new(1)).ok();
}
