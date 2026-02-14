use std::fmt::Write;

use bezier_nd::{bignum::RationalN, Num};

pub fn matrix_to_strings<F: std::fmt::Display>(m: &[F], num_cols: usize) -> Vec<String> {
    let mut result = vec![];
    for row in m.chunks_exact(num_cols) {
        let mut line: String = "|".into();
        for p in row {
            line.push_str(&format!(" {p:16}"));
        }
        line.push('|');
        result.push(line);
    }
    result
}

pub fn eprintln_matrix<F: std::fmt::Display>(m: &[F], num_cols: usize) {
    let r = matrix_to_strings(m, num_cols);
    for r in r {
        eprintln!("{r}");
    }
}

pub fn eprintln_rational_matrix<const N: usize>(m: &[RationalN<N>], num_cols: usize) {
    let (i, lcm) = RationalN::<_>::with_common_denom(m.iter());
    let m2: Vec<_> = i.collect();
    eprintln!("1/{lcm}*");
    eprintln_matrix(&m2, num_cols);
}

pub fn eprintln_f64_matrix<const N: usize>(m: &[RationalN<N>], num_cols: usize) {
    let (i, lcm) = RationalN::<_>::with_common_denom(m.iter());
    let m2: Vec<_> = i
        .map(|i| {
            let f: RationalN<N> = (i, lcm).into();
            f64::try_from(&f).unwrap()
        })
        .collect();

    eprintln_matrix(&m2, num_cols);
}

pub fn string_of_f64_matrix<const N: usize>(m: &[RationalN<N>]) -> String {
    let (i, lcm) = RationalN::<_>::with_common_denom(m.iter());

    i.fold(format!("["), |mut s, i| {
        let f: RationalN<N> = (i, lcm).into();
        let f = f64::try_from(&f).unwrap();
        s.push_str(&format!(" {f:?},"));
        s
    }) + "]"
}

pub fn string_of_f64_matrix_f64(m: &[f64]) -> String {
    m.iter().fold(format!("["), |mut s, f| {
        s.push_str(&format!(" {f:?},"));
        s
    }) + "]"
}
