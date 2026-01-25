//a Imports
mod utils;
use bezier_nd::bernstein::bezier_fns::generate_elevate_by_one_matrix;
use bezier_nd::{Bezier, BezierND, Float};
use geo_nd::{
    matrix,
    vector::{self},
    FArray,
};
use utils::test_beziers_approx_eq;
use utils::{assert_near_equal, assert_near_identity, float_iter, generate_bernstein_matrix};

#[test]
fn bisect() {
    let p0: FArray<f64, 2> = [0., 0.].into();
    let p1: FArray<f64, 2> = [10., 0.].into();
    let p2: FArray<f64, 2> = [6., 1.].into();
    let p3: FArray<f64, 2> = [20., 5.].into();

    test_beziers_approx_eq(
        &Bezier::cubic(&p0, &p1, &p2, &p3),
        &Bezier::cubic(&p0, &p1, &p2, &p3),
    );
    test_beziers_approx_eq(
        &Bezier::cubic(&p2, &p1, &p3, &p0),
        &Bezier::cubic(&p2, &p1, &p3, &p0),
    );
    test_beziers_approx_eq(
        &Bezier::quadratic(&p2, &p1, &p3),
        &Bezier::quadratic(&p2, &p1, &p3),
    );
    test_beziers_approx_eq(
        &Bezier::quadratic(&p0, &p1, &p2),
        &Bezier::quadratic(&p0, &p1, &p2),
    );
    test_beziers_approx_eq(&Bezier::line(&p0, &p2), &Bezier::line(&p0, &p2));
}

#[test]
fn bernstein_matrix() {
    let mut bern_n = [0.0f64; 100];
    generate_bernstein_matrix(&mut bern_n[0..3], 2, &[0.]);
    assert_eq!(&bern_n[0..3], &[1., 0., 0.]);
    generate_bernstein_matrix(&mut bern_n[0..3], 2, &[0.5]);
    assert_eq!(&bern_n[0..3], &[0.25, 0.5, 0.25]);
    generate_bernstein_matrix(&mut bern_n[0..3], 2, &[1.0]);
    assert_eq!(&bern_n[0..3], &[0., 0., 1.]);
}

fn generate_bs_reduce_matrix(
    degree: usize,
    ts: &[f64],
    reduce: &mut [f64],
    elevated_reduce: &mut [f64],
) {
    assert_eq!(
        ts.len(),
        degree + 1,
        "Generated a reduction matrix to {degree} requires {} points",
        degree + 1
    );
    let n2 = (degree + 1) * (degree + 1);
    // bern_n is a (degree+1) * (degree+1) matrix
    let mut bern_n = [0.0f64; 100];
    let mut bern_np1 = [0.0f64; 100];
    let mut tr0 = [0.0f64; 100];
    let mut tr1 = [0.0f64; 100];
    let mut pivot = [0; 100];
    let mut lu = [0.0f64; 100];
    generate_bernstein_matrix(&mut bern_n[0..n2], degree, &ts);
    generate_bernstein_matrix(
        &mut bern_np1[0..(degree + 2) * (degree + 1)],
        degree + 1,
        &ts,
    );
    let mut bern_n_inverse = bern_n.clone();
    assert_ne!(
        matrix::lup_decompose(degree + 1, &bern_n[0..n2], &mut lu[0..n2], &mut pivot),
        0.0,
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
    let mut test = [0.0f64; 100];
    let mut elevate = [0.0f64; 100];
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

    eprintln!("Bernstein reduce for ts {ts:?} degree {degree}");
    eprintln!("  &{:?}", &reduce[0..(degree + 1) * (degree + 2)]);
    eprintln!("Bernstein elevated reduce for ts {ts:?} degree {degree}");
    eprintln!("  &{:?}", &elevated_reduce[0..(degree + 2) * (degree + 2)]);
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
    let mut elevate = [0.0f64; 100];
    let scale = generate_elevate_by_one_matrix(&mut elevate, 2);
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

    generate_bs_reduce_matrix(1, &[0., 1.0], &mut reduce, &mut elevated_reduce);
    assert_near_equal(&reduce, &bezier_nd::REDUCE_BY_ONE_BS_UNIFORM_F64[0]);
    assert_near_equal(
        &elevated_reduce,
        &bezier_nd::ELEVATED_REDUCE_BY_ONE_BS_UNIFORM_F64[0],
    );

    generate_bs_reduce_matrix(2, &[0., 0.5, 1.0], &mut reduce, &mut elevated_reduce);
    assert_near_equal(&reduce, &bezier_nd::REDUCE_BY_ONE_BS_UNIFORM_F64[1]);
    assert_near_equal(
        &elevated_reduce,
        &bezier_nd::ELEVATED_REDUCE_BY_ONE_BS_UNIFORM_F64[1],
    );

    generate_bs_reduce_matrix(
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
        5,
        &[0., 1.0 / 5.0, 2.0 / 5.0, 3.0 / 5.0, 4.0 / 5.0, 1.0],
        &mut reduce,
        &mut elevated_reduce,
    );
    assert_near_equal(&reduce, &bezier_nd::REDUCE_BY_ONE_BS_UNIFORM_F64[4]);

    generate_bs_reduce_matrix(
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
    let dm = b.metric_dm_est(&nb, 1000);
    let df = b.metric_df(&nb);
    let dc = b.metric_dc(&nb);
    eprintln!("Metrics for b/nb: dm {dm} < dc {dc} < df {df}");
    assert!(dm < dc);
    assert!(dc < df);

    let (b0, b1) = b.bisect();
    for t in float_iter(-1.0, 1.0, 100) {
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
    let dm0 = b0.metric_dm_est(&nb0, 1000);
    let df0 = b0.metric_df(&nb0);
    let dc0 = b0.metric_dc(&nb0);
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
        for tp in float_iter(0., 1., 100) {
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

use bezier_nd::BezierBuilder;
fn build_bezier<F: Float, const N: usize, const D: usize>(
    builder: BezierBuilder<F, D>,
) -> BezierND<F, N, D> {
    let degree = builder.bezier_min_degree();
    let n2 = (degree + 1) * (degree + 1);
    let mut bern_n = [F::zero(); 100];

    let mut basis = vec![];
    let mut pts = vec![];
    for c in builder.iter() {
        let t = c.at();
        let pt = c.posn();
        pts.push(*pt);
        generate_bernstein_matrix(&mut bern_n[0..(degree + 1)], degree, &[t]);
        basis.extend(bern_n.iter().take(degree + 1));
    }

    assert_eq!(
        pts.len(),
        degree + 1,
        "Degree must match the number of constraints given"
    );
    let bezier = BezierND::<F, N, D>::new(&pts);

    let mut basis_inverse = basis.clone();
    let mut lu = basis.clone();
    let mut pivot = vec![0; degree + 1];
    assert_ne!(
        matrix::lup_decompose(degree + 1, &basis[0..n2], &mut lu[0..n2], &mut pivot),
        F::zero(),
        "Matrix is invertible"
    );

    let mut tr0 = vec![F::zero(); degree + 1];
    let mut tr1 = vec![F::zero(); degree + 1];
    assert!(
        matrix::lup_invert(
            degree + 1,
            &lu,
            &pivot,
            &mut basis_inverse,
            &mut tr0,
            &mut tr1
        ),
        "Matrix is invertible"
    );
    bezier.apply_matrix(&basis_inverse, degree)
}

#[test]
fn build() {
    let mut builder = BezierBuilder::default();
    builder.add_point_at(0.0_f32, [1., 2.]);
    builder.add_point_at(0.5_f32, [3., 0.]);

    let b = build_bezier::<_, 10, 2>(builder);
    dbg!(b);
    //    assert!(false);
}
