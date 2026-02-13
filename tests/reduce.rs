mod utils;
use bezier_nd::bernstein_fns::elevate_reduce_matrix::generate_elevate_by_one_matrix;
use bezier_nd::bignum::RationalN;
use bezier_nd::{metrics, BezierEval, BezierND, BezierReduce, BezierReduction, BezierSplit, Num};
use geo_nd::{matrix, vector};
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

/// Generate a Bernstein reduction matrix *to* degree (using ts provided)
///
/// ts must be of length degree+1, and the reduction will keep the values
/// of the reduced Bezier at these values of t *unchanged* compared to the
/// original Bezier
fn generate_bs_reduce_matrix<N: Num + Into<f64>>(
    reduce_from_degree: usize,
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
    let mut bern_n = vec![N::ZERO; n2];
    let mut bern_np1 = vec![N::ZERO; (degree + 1) * (reduce_from_degree + 1)];
    let mut tr0 = vec![N::ZERO; degree + 1];
    let mut tr1 = vec![N::ZERO; degree + 1];
    let mut pivot = vec![0; degree + 1];
    let mut lu = bern_n.clone();

    // bern_n is the matrix that can be applied to a Bezier of 'degree' to
    // generate the points at the values of t given in 'ts'
    //
    // bern_n is (degree+1) * (degree+1)
    generate_bernstein_matrix_br(&mut bern_n[0..n2], degree, &ts);

    // bern_np1 is the matrix that can be applied to a Bezier of 'reduce_from_degree' to
    // generate the points at the values of t given in 'ts'
    //
    // bern_np1 is (degree+1) * (reduce_from_degree+1)
    generate_bernstein_matrix_br(
        &mut bern_np1[0..(reduce_from_degree + 1) * (degree + 1)],
        reduce_from_degree,
        &ts,
    );

    // bern_n_inverse is such that bern_n_inverse * bern_n == Identity
    //
    // bern_n_inverse is (degree+1) * (degree+1)
    let mut bern_n_inverse = bern_n.clone();
    assert_ne!(
        matrix::lup_decompose(degree + 1, &bern_n, &mut lu, &mut pivot),
        N::ZERO,
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

    // reduce = bern_n_inverse * bern_np1
    //
    // (degree+1) * (degree+1) times (degree+1) * (reduce_from_degree+1) yields...
    //
    // reduce is (degree+1) *(reduce_from_degree+1)
    matrix::multiply_dyn(
        degree + 1,
        degree + 1,
        reduce_from_degree + 1,
        &bern_n_inverse,
        &bern_np1,
        reduce,
    );

    // Find reduction of elevation
    let mut test = [N::zero(); 100];
    let elevate = generate_elevate_by_n_matrix(degree, reduce_from_degree - degree);

    // Validate that reduce of elevation is the identity (approx)
    matrix::multiply_dyn(
        degree + 1,
        reduce_from_degree + 1,
        degree + 1,
        &reduce,
        &elevate,
        &mut test,
    );
    eprintln!(
        "test {:?}",
        test.iter()
            .take((degree + 1) * (degree + 1))
            .map(|s| s.to_string())
            .collect::<Vec<_>>()
    );
    assert_near_identity(degree + 1, &test);

    // Generate the elevated_reduce matrix
    //
    matrix::multiply_dyn(
        reduce_from_degree + 1,
        degree + 1,
        reduce_from_degree + 1,
        &elevate,
        &reduce,
        elevated_reduce,
    );

    eprintln!("Bernstein reduce for ts from {reduce_from_degree} to degree {degree}");
    eprintln!(
        "  &[{}]",
        reduce
            .iter()
            .take((degree + 1) * (reduce_from_degree + 1))
            .fold(String::new(), |string, s| string
                + &(<N as Into<f64>>::into(*s)).to_string()
                + ", ")
    );
    eprintln!(
        "  &[{}]",
        reduce
            .iter()
            .take((degree + 1) * (reduce_from_degree + 1))
            .fold(String::new(), |string, s| string + &s.to_string() + ", ")
    );
    eprintln!("Bernstein elevated reduce for ts {reduce_from_degree} to degree {degree}");
    eprintln!(
        "  &[{}]",
        elevated_reduce
            .iter()
            .take((reduce_from_degree + 1) * (reduce_from_degree + 1))
            .fold(String::new(), |string, s| string
                + &(<N as Into<f64>>::into(*s)).to_string()
                + ", ")
    );
    eprintln!(
        "  &[{}]",
        elevated_reduce
            .iter()
            .take((reduce_from_degree + 1) * (reduce_from_degree + 1))
            .fold(String::new(), |string, s| string + &s.to_string() + ", ")
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
        generate_bs_reduce_matrix(
            (degree + 1) as usize,
            degree as usize,
            &ts,
            &mut reduce,
            &mut elevated_reduce,
        );
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
}

#[test]
fn validate_bs_reductions_to_quad() {
    let mut reduce = [RationalN::<8>::default(); 200];
    let mut elevated_reduce = [RationalN::<8>::default(); 200];

    for degree in [1, 2, 3] {
        let ts: Vec<RationalN<8>> = (0..(degree + 1))
            .map(|i| (i as i64, degree).into())
            .collect();

        for reduce_from_degree in (degree + 1)..11 {
            generate_bs_reduce_matrix(
                reduce_from_degree as usize,
                degree as usize,
                &ts,
                &mut reduce,
                &mut elevated_reduce,
            );
            let (i, _lcm) = RationalN::<8>::with_common_denom(reduce.iter());
            eprintln!(
                "\n\n****\nreduce = [{}]",
                i.fold(String::new(), |s, a| s + &a.to_string() + ", ")
            );
            let (i, lcm) = RationalN::<8>::with_common_denom(elevated_reduce.iter());
            eprintln!(
                "elevated_reduce = {lcm} : [{}]",
                i.fold(String::new(), |s, a| s + &a.to_string() + ", ")
            );
        }
    }
    // assert!(false);
}

#[test]
fn bernstein_reduce_matrix_c_to_q() {
    let mut bern_n = [0.0f64; 100];
    let mut bern_np1 = [0.0f64; 100];

    // Find reduction from cubic (degree 3, 4 pts) to quadratic using ts of 0, 1/2, 1
    let ts = [0.0, 0.5, 1.0];
    generate_bernstein_matrix(&mut bern_n[0..9], 2, &ts);
    assert_near_equal(&bern_n, &[1.0, 0., 0., 0.25, 0.5, 0.25, 0., 0., 1.]);
    generate_bernstein_matrix(&mut bern_np1[0..12], 3, &ts);
    assert_near_equal(
        &bern_np1,
        &[1.0, 0., 0., 0., 0.125, 0.375, 0.375, 0.125, 0., 0., 0., 1.],
    );
    let mut bern_n_inverse = bern_n.clone();
    assert!(
        matrix::invert::<_, 9, 3>(&mut bern_n_inverse.first_chunk_mut::<9>().unwrap()),
        "Bernstein basis must be invertible"
    );
    assert_near_equal(&bern_n_inverse, &[1.0, 0., 0., -0.5, 2.0, -0.5, 0., 0., 1.]);

    // Reduction matrix is 3x4
    let reduce = matrix::multiply::<_, 9, 12, 12, 3, 3, 4>(
        bern_n_inverse.first_chunk::<9>().unwrap(),
        bern_np1.first_chunk::<12>().unwrap(),
    );
    eprintln!("bernstein_reduce[3] = {reduce:?}");
    assert_near_equal(
        &reduce,
        &[
            1.0, 0.0, 0., 0., -0.25, 0.75, 0.75, -0.25, 0.0, 0.0, 0., 1.0,
        ],
    );

    // Find elevation from quadratic to cubic
    let (scale, mut elevate) = generate_elevate_by_one_matrix(2);
    for e in elevate.iter_mut() {
        *e /= scale;
    }

    // Find reduction of elevation
    let reduced_elevate = matrix::multiply::<_, 12, 12, 9, 3, 4, 3>(
        reduce.first_chunk::<12>().unwrap(),
        elevate.first_chunk::<12>().unwrap(),
    );
    dbg!(&reduced_elevate);
    assert_near_equal(&reduced_elevate, &matrix::identity3());
}

#[test]
fn bernstein_reduce_matrix() {
    let mut reduce = [0.0; 100];
    let mut elevated_reduce = [0.0; 100];

    generate_bs_reduce_matrix(2, 1, &[0., 1.0], &mut reduce, &mut elevated_reduce);
    assert_near_equal(&reduce, &bezier_nd::REDUCE_BY_ONE_BS_UNIFORM_F64[0]);
    assert_near_equal(
        &elevated_reduce,
        &bezier_nd::ELEVATED_REDUCE_BY_ONE_BS_UNIFORM_F64[0],
    );

    generate_bs_reduce_matrix(3, 2, &[0., 0.5, 1.0], &mut reduce, &mut elevated_reduce);
    assert_near_equal(&reduce, &bezier_nd::REDUCE_BY_ONE_BS_UNIFORM_F64[1]);
    assert_near_equal(
        &elevated_reduce,
        &bezier_nd::ELEVATED_REDUCE_BY_ONE_BS_UNIFORM_F64[1],
    );

    generate_bs_reduce_matrix(
        4,
        3,
        &[0., 1.0 / 3.0, 2.0 / 3.0, 1.0],
        &mut reduce,
        &mut elevated_reduce,
    );
    assert_near_equal(&reduce, &bezier_nd::REDUCE_BY_ONE_BS_UNIFORM_F64[2]);
    assert_near_equal(
        &elevated_reduce,
        &bezier_nd::ELEVATED_REDUCE_BY_ONE_BS_UNIFORM_F64[2],
    );

    generate_bs_reduce_matrix(
        5,
        4,
        &[0., 1.0 / 4.0, 2.0 / 4.0, 3.0 / 4.0, 1.0],
        &mut reduce,
        &mut elevated_reduce,
    );
    assert_near_equal(&reduce, &bezier_nd::REDUCE_BY_ONE_BS_UNIFORM_F64[3]);
    assert_near_equal(
        &elevated_reduce,
        &bezier_nd::ELEVATED_REDUCE_BY_ONE_BS_UNIFORM_F64[3],
    );

    generate_bs_reduce_matrix(
        6,
        5,
        &[0., 1.0 / 5.0, 2.0 / 5.0, 3.0 / 5.0, 4.0 / 5.0, 1.0],
        &mut reduce,
        &mut elevated_reduce,
    );
    assert_near_equal(&reduce, &bezier_nd::REDUCE_BY_ONE_BS_UNIFORM_F64[4]);

    generate_bs_reduce_matrix(
        7,
        6,
        &[
            0.,
            1.0 / 6.0,
            2.0 / 6.0,
            3.0 / 6.0,
            4.0 / 6.0,
            5.0 / 6.0,
            1.0,
        ],
        &mut reduce,
        &mut elevated_reduce,
    );
}

#[test]
fn reduce_and_elevate_cubic_in_parts() {
    let p0 = [0., 0.];
    let p1 = [10., 0.];
    let p2 = [6., 1.];
    let p3 = [20., 5.];

    let b = BezierND::<f64, 4, 2>::new(&[p0, p1, p2, p3]);
    let nb = b.apply_matrix(&bezier_nd::ELEVATED_REDUCE_BY_ONE_BS_UNIFORM_F64[1], 3);
    let ts = [0., 0.5, 1.0];
    for t in ts {
        assert!(
            vector::distance_sq(&b.point_at(t), &nb.point_at(t)) < 1E-4,
            "Uniform reduction to quadratic has same point for t {t} in {ts:?}"
        );
    }
    let dm = metrics::dm_sq_est(b.control_points(), nb.control_points(), 1000)
        .unwrap()
        .sqrt();
    let df = metrics::df_sq(b.control_points(), nb.control_points())
        .unwrap()
        .sqrt();
    let dc = metrics::dc_sq(b.control_points(), nb.control_points())
        .unwrap()
        .sqrt();
    eprintln!("Metrics for b/nb: dm {dm} < dc {dc} < df {df}");
    assert!(dm < dc);
    assert!(dc < df);

    let (b0, b1) = b.split();
    for t in float_iter_between(-1.0, 1.0, 100) {
        let p = b.point_at(t);
        let p0 = b0.point_at(t * 2.0);
        let p1 = b1.point_at(t * 2.0 - 1.0);
        assert!(
            vector::distance_sq(&p, &p0) < 1E-4,
            "Bisection first half should match Bezier with t*=2 at {t} {p:?} {p0:?}"
        );
        assert!(
            vector::distance_sq(&p, &p1) < 1E-4,
            "Bisection second half should match Bezier with t=2t-1 at {t} {p:?} {p0:?}"
        );
    }

    let nb0 = b0.apply_matrix(&bezier_nd::ELEVATED_REDUCE_BY_ONE_BS_UNIFORM_F64[1], 3);
    let dm0 = metrics::dm_sq_est(b0.control_points(), nb0.control_points(), 1000)
        .unwrap()
        .sqrt();
    let df0 = metrics::df_sq(b0.control_points(), nb0.control_points())
        .unwrap()
        .sqrt();
    let dc0 = metrics::dc_sq(b0.control_points(), nb0.control_points())
        .unwrap()
        .sqrt();
    eprintln!("Metrics for b0/nb0: dm {dm0} < dc {dc0} < df {df0}");
    assert!(dm0 < dc0);
    assert!(dc0 < df0);

    assert!(dm0 < dm);
    assert!(dc0 < dc);
    assert!(df0 < df);
}

#[test]
fn reduce_and_elevate_cubic() {
    let p0 = [0., 0.];
    let p1 = [10., -1.];
    let p2 = [6., 1.];
    let p3 = [20., 5.];
    let max_dc = 0.03;

    let b = BezierND::<f64, 4, 2>::new(&[p0, p1, p2, p3]);
    assert_eq!(b.degree(), 3);
    let b_split = b.reduce_and_split_iter(
        &bezier_nd::REDUCE_BY_ONE_BS_UNIFORM_F64[1],
        &bezier_nd::ELEVATED_REDUCE_BY_ONE_BS_UNIFORM_F64[1],
        2,
        max_dc,
    );
    let mut t = 0.0;
    for (ns, bs) in b_split {
        let t0 = t;
        let dt = 1.0 / (2.0_f64).powi(ns as i32);
        let t1 = t + dt;
        let p0 = b.point_at(t0);
        let p1 = b.point_at(t1);
        let ps0 = bs.point_at(0.0);
        let ps1 = bs.point_at(1.0);
        assert!(
            vector::distance_sq(&p0, &ps0) < 1E-4,
            "Split should match at pt0 {t0} {p0:?} {ps0:?}"
        );
        assert!(
            vector::distance_sq(&p1, &ps1) < 1E-4,
            "Split should match at pt1 {t1} {p1:?} {ps1:?}"
        );
        for tp in float_iter(100) {
            let p = b.point_at(t0 + tp * dt);
            let ps = bs.point_at(tp);
            assert!(
                vector::distance_sq(&p, &ps) < max_dc,
                "Split should be no more than {max_dc} apartat pt {tp} of {t0}->{t1} {p:?} {ps:?}"
            );
        }
        t = t1;
    }
}

fn test_reduce_least_squares<F, const D: usize, B>(bezier: &B)
where
    B: BezierReduce<F, [F; D]>,
    F: Num,
{
    assert!(
        bezier.can_reduce(BezierReduction::LeastSquares),
        "Must be able to reduce by least squares"
    );
}

#[test]
fn least_squares() {
    test_reduce_least_squares(&[[1.0_f32], [1.0_f32], [1.0_f32]]);
}
