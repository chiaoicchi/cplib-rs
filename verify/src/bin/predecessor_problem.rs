use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::collections::int_set::IntSet;

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
    let t: &[u8] = iter.next().unwrap();
    let mut int_set = IntSet::new(n);
    for (i, &t) in t.iter().enumerate() {
        if t == b'1' {
            int_set.insert(i);
        }
    }
    for _ in 0..q {
        let c = parse!(u8);
        let k = parse!(usize);
        match c {
            0 => {
                int_set.insert(k);
            }
            1 => {
                int_set.remove(k);
            }
            2 => {
                let ans = int_set.contains(k);
                writeln!(stdout, "{}", ans as u8).ok();
            }
            3 => {
                let ans = int_set.floor(k).unwrap_or(!0);
                writeln!(stdout, "{}", ans as isize).ok();
            }
            4 => {
                let ans = int_set.ceil(k).unwrap_or(!0);
                writeln!(stdout, "{}", ans as isize).ok();
            }
            _ => unreachable!(),
        }
    }
}
