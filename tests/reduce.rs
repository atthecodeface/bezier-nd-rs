mod utils;
use bezier_nd::bernstein::bezier_fns::generate_elevate_by_one_matrix;
use geo_nd::matrix;
use num::rational::Rational64;
use utils::*;

#[test]
fn bernstein_matrix_br() {
    let mut bern_n = [Rational64::default(); 100];
    generate_bernstein_matrix_br(&mut bern_n[0..3], 2, &[0_i64.into()]);
    assert_eq!(bern_n[0], 1_i64.into());
    assert_eq!(bern_n[1], 0_i64.into());
    assert_eq!(bern_n[2], 0_i64.into());
    generate_bernstein_matrix_br(&mut bern_n[0..3], 2, &[(1_i64, 2_i64).into()]);
    assert_eq!(bern_n[0], (1_i64, 4_i64).into());
    assert_eq!(bern_n[1], (1_i64, 2_i64).into());
    assert_eq!(bern_n[2], (1_i64, 4_i64).into());
    generate_bernstein_matrix_br(&mut bern_n[0..3], 2, &[1_i64.into()]);
    assert_eq!(bern_n[0], 0_i64.into());
    assert_eq!(bern_n[1], 0_i64.into());
    assert_eq!(bern_n[2], 1_i64.into());
}

fn generate_bs_reduce_matrix(
    degree: usize,
    ts: &[Rational64],
    reduce: &mut [Rational64],
    elevated_reduce: &mut [Rational64],
) {
    assert_eq!(
        ts.len(),
        degree + 1,
        "Generated a reduction matrix to {degree} requires {} points",
        degree + 1
    );
    let n2 = (degree + 1) * (degree + 1);
    // bern_n is a (degree+1) * (degree+1) matrix
    let mut bern_n = [Rational64::default(); 200];
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
    let mut test = [Rational64::default(); 100];
    let mut elevate = [Rational64::default(); 100];
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

    /*
    matrix::multiply_dyn(
        degree + 2,
        degree + 1,
        degree + 2,
        &elevate,
        &reduce,
        elevated_reduce,
    );

    eprintln!("Bernstein reduce for ts {ts:?} degree {degree}");
    eprintln!("  &{:?}", &reduce[0..(degree + 1) * (degree + 2)]);
    eprintln!("Bernstein elevated reduce for ts {ts:?} degree {degree}");
    eprintln!("  &{:?}", &elevated_reduce[0..(degree + 2) * (degree + 2)]);
    */
}
#[test]
fn bernstein_reduce_matrix() {
    let mut reduce = [Rational64::default(); 200];
    let mut elevated_reduce = [Rational64::default(); 200];

    let degree: i64 = 6;
    let ts: Vec<Rational64> = (0..(degree + 1)).map(|i: i64| (i, degree).into()).collect();
    generate_bs_reduce_matrix(degree as usize, &ts, &mut reduce, &mut elevated_reduce);
    //   assert_near_equal(&reduce, &bezier_nd::REDUCE_BY_ONE_BS_UNIFORM_F64[0]);
    //   assert_near_equal(
    //       &elevated_reduce,
    //      &bezier_nd::ELEVATED_REDUCE_BY_ONE_BS_UNIFORM_F64[0],
    //  );
}
