mod utils;
use bezier_nd::bernstein::bezier_fns::generate_elevate_by_one_matrix;
use bezier_nd::bignum::RationalN;
use geo_nd::{matrix, Num};
use utils::*;

#[test]
fn bernstein_matrix_br() {
    type N = RationalN<8>;
    let mut bern_n = [N::default(); 100];
    generate_bernstein_matrix_br(&mut bern_n[0..3], 2, &[0_i64.into()]);
    assert_eq!(bern_n[0], 1_i64.into());
    assert_eq!(bern_n[1], 0_i64.into());
    assert_eq!(bern_n[2], 0_i64.into());
    generate_bernstein_matrix_br(&mut bern_n[0..3], 2, &[(1_i64, 2_u64).into()]);
    assert_eq!(bern_n[0], (1_i64, 4_u64).into());
    assert_eq!(bern_n[1], (1_i64, 2_u64).into());
    assert_eq!(bern_n[2], (1_i64, 4_u64).into());
    generate_bernstein_matrix_br(&mut bern_n[0..3], 2, &[1_i64.into()]);
    assert_eq!(bern_n[0], 0_i64.into());
    assert_eq!(bern_n[1], 0_i64.into());
    assert_eq!(bern_n[2], 1_i64.into());
}

fn generate_bs_reduce_matrix<N: Num + From<i64> + Into<f64>>(
    degree: usize,
    ts: &[N],
    reduce: &mut [N],
    elevated_reduce: &mut [N],
) {
    assert_eq!(
        ts.len(),
        degree + 1,
        "Generated a reduction matrix to {degree} requires {} points",
        degree + 1
    );
    let n2 = (degree + 1) * (degree + 1);
    // bern_n is a (degree+1) * (degree+1) matrix
    let mut bern_n = [N::ZERO; 200];
    let mut bern_np1 = bern_n.clone();
    let mut tr0 = bern_n.clone();
    let mut tr1 = bern_n.clone();
    let mut pivot = [0; 200];
    let mut lu = bern_n.clone();
    generate_bernstein_matrix_br(&mut bern_n[0..n2], degree, &ts);
    generate_bernstein_matrix_br(
        &mut bern_np1[0..(degree + 2) * (degree + 1)],
        degree + 1,
        &ts,
    );
    let mut bern_n_inverse = bern_n.clone();

    assert_ne!(
        matrix::lup_decompose(degree + 1, &bern_n[0..n2], &mut lu[0..n2], &mut pivot),
        (0_i64).into(),
        "Matrix is invertible"
    );
    assert!(
        matrix::lup_invert(
            degree + 1,
            &lu,
            &pivot,
            &mut bern_n_inverse,
            &mut tr0,
            &mut tr1
        ),
        "Matrix is invertible"
    );
    // Reduction matrix is 3x4
    matrix::multiply_dyn(
        degree + 1,
        degree + 1,
        degree + 2,
        &bern_n_inverse,
        &bern_np1,
        reduce,
    );

    // Find reduction of elevation
    let mut test = [N::zero(); 100];
    let mut elevate = [N::zero(); 100];
    let scale = generate_elevate_by_one_matrix(&mut elevate, degree);
    for e in elevate.iter_mut() {
        *e /= scale;
    }
    matrix::multiply_dyn(
        degree + 1,
        degree + 2,
        degree + 1,
        &reduce,
        &elevate,
        &mut test,
    );
    assert_near_identity(degree + 1, &test);

    matrix::multiply_dyn(
        degree + 2,
        degree + 1,
        degree + 2,
        &elevate,
        &reduce,
        elevated_reduce,
    );

    eprintln!("Bernstein reduce for ts degree {degree}");
    eprintln!(
        "  &[{}]",
        reduce
            .iter()
            .take((degree + 1) * (degree + 2))
            .fold(String::new(), |string, s| string
                + &(<N as Into<f64>>::into(*s)).to_string()
                + ", ")
    );
    eprintln!("Bernstein elevated reduce for ts degree {degree}");
    eprintln!(
        "  &[{}]",
        elevated_reduce
            .iter()
            .take((degree + 2) * (degree + 2))
            .fold(String::new(), |string, s| string
                + &(<N as Into<f64>>::into(*s)).to_string()
                + ", ")
    );
}

#[test]
fn validate_bernstein_reduce_matrices() {
    let mut reduce = [RationalN::<8>::default(); 200];
    let mut elevated_reduce = [RationalN::<8>::default(); 200];

    for degree in 1..7 {
        let ts: Vec<RationalN<8>> = (0..(degree + 1))
            .map(|i| (i as i64, degree).into())
            .collect();
        generate_bs_reduce_matrix(degree as usize, &ts, &mut reduce, &mut elevated_reduce);
        let reduce_f64: Vec<f64> = reduce.iter().map(|v| v.into()).collect();
        assert_near_equal(
            &bezier_nd::REDUCE_BY_ONE_BS_UNIFORM_F64[degree as usize - 1],
            &reduce_f64,
        );
        let elevated_reduce_f64: Vec<f64> = elevated_reduce.iter().map(|v| v.into()).collect();
        assert_near_equal(
            &bezier_nd::ELEVATED_REDUCE_BY_ONE_BS_UNIFORM_F64[degree as usize - 1],
            &elevated_reduce_f64,
        );
    }
    // assert!(false);
}
