use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::algebra::canonical::{Affine, Canonical};
use cplib::collections::foldable_queue::FoldableQueue;
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

    let q = parse!(u32);
    let mut queue = FoldableQueue::new(Affine(Canonical::new()));

    for _ in 0..q {
        let t = parse!(u8);
        match t {
            0 => {
                let a = parse!(u32);
                let b = parse!(u32);
                queue.push((fp!(a, mod P), fp!(b)));
            }
            1 => {
                queue.pop();
            }
            2 => {
                let x = parse!(u32);
                let ans = match queue.fold() {
                    Some((a, b)) => a * fp!(x) + b,
                    None => fp!(x),
                };
                writeln!(stdout, "{ans}").ok();
            }
            _ => unreachable!(),
        }
    }
}
