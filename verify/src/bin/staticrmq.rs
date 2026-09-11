use std::io::{BufWriter, Read, Write, stdin, stdout};

use cplib::algebra::min::Min;
use cplib::collections::sparse_table::SparseTable;

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

    let a: Vec<u32> = (0..n).map(|_| parse!(u32)).collect();
    let sparse_table = SparseTable::from_vec(Min::new(), a);

    for _ in 0..q {
        let l = parse!(usize);
        let r = parse!(usize);
        writeln!(stdout, "{}", sparse_table.fold(l..r).unwrap()).ok();
    }
}
