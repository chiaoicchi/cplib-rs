use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::collections::range_set::RangeSet;

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
    let mut range_set = RangeSet::new();
    for (i, &t) in t.iter().enumerate() {
        if t == b'1' {
            range_set.insert(i, i + 1);
        }
    }

    for _ in 0..q {
        let c = parse!(u8);
        let k = parse!(usize);
        match c {
            0 => range_set.insert(k, k + 1),
            1 => range_set.remove(k, k + 1),
            2 => {
                let ans = range_set.contains(k);
                writeln!(stdout, "{}", ans as u8).ok();
            }
            3 => {
                let ans = range_set
                    .range(k, n)
                    .next()
                    .map_or(-1, |(a, _)| a.max(k) as i64);
                writeln!(stdout, "{ans}").ok();
            }
            4 => {
                let ans = range_set
                    .range(0, k + 1)
                    .next_back()
                    .map_or(-1, |(_, b)| (b - 1).min(k) as i64);
                writeln!(stdout, "{ans}").ok();
            }
            _ => unreachable!(),
        }
    }
}
