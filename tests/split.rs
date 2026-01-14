//a Imports
use bezier_nd::Bezier;
use geo_nd::{
    matrix,
    vector::{self, reduce},
    FArray,
};

fn float_iter(t0: f64, t1: f64, n: usize) -> impl Iterator<Item = f64> {
    assert!(
        n >= 2,
        "Float iterator must have at least two steps for begin and end"
    );
    let r = t1 - t0;
    assert!(r >= 0.0, "Float range must be positive");
    (0..n).map(move |i: usize| (i as f64) / (r * ((n - 1) as f64)))
}

fn test_subsection(b: &Bezier<f64, 2>, sub: &Bezier<f64, 2>, t0: f64, t1: f64) {
    eprintln!("Testing subsections of beziers {b} {sub} {t0} {t1}");
    for sub_t in float_iter(0.0, 1.0, 100) {
        let t = t0 + (t1 - t0) * sub_t;
        let p = b.point_at(t);
        let sub_p = sub.point_at(sub_t);
        let d2 = vector::distance_sq(&p, &sub_p);
        assert!(d2<1E-4,
        "Points at bezier {t} : {p:?} and subbezier {sub_t} : {sub_p:?} should be roughly the same but have distance {d2}");
    }
}

fn test_beziers_approx_eq(b0: &Bezier<f64, 2>, b1: &Bezier<f64, 2>) {
    test_subsection(b0, b1, 0.0, 1.0);
    test_subsection(b0, &b1.bezier_between(0.0, 1.0), 0.0, 1.0);
    test_subsection(b0, &b1.bezier_between(0.1, 0.4), 0.1, 0.4);
    test_subsection(b0, &b1.bisect().0, 0.0, 0.5);
    test_subsection(b0, &b1.bisect().1, 0.5, 1.0);
}

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
fn elevate() {
    let p0: FArray<f64, 2> = [0., 0.].into();
    let p1: FArray<f64, 2> = [10., 0.].into();
    let p2: FArray<f64, 2> = [6., 1.].into();
    let p3: FArray<f64, 2> = [20., 5.].into();

    let b = Bezier::quadratic(&p0, &p1, &p2);
    let b2 = b.clone().elevate();
    test_beziers_approx_eq(&b, &b2);

    let b = Bezier::quadratic(&p3, &p1, &p2);
    let b2 = b.clone().elevate();
    test_beziers_approx_eq(&b, &b2);

    let b = Bezier::line(&p3, &p1);
    let b2 = b.clone().elevate();
    test_beziers_approx_eq(&b, &b2);
}

fn generate_elevate_by_one_matrix(matrix: &mut [f64], degree: usize) {
    assert!(
        matrix.len() >= (degree + 1) * (degree + 2),
        "Must have enough room in matrix for coeffs for degree -> degreee+1"
    );
    for i in 0..(degree + 1) {
        for j in 0..(degree + 2) {
            matrix[j * (degree + 1) + i] =
                bezier_nd::bezier::Bezier::<_, 2, 10>::elevation_by_one_matrix_ele(degree, i, j);
        }
    }
}

fn generate_bernstein_matrix(matrix: &mut [f64], degree: usize, ts: &[f64]) {
    assert_eq!(
        matrix.len(),
        (degree + 1) * ts.len(),
        "Trying to generate an {} by {} matrix, so must be given a matrix of that size",
        ts.len(),
        degree + 1
    );
    for (r, t) in ts.iter().enumerate() {
        for c in 0..(degree + 1) {
            matrix[r * (degree + 1) + c] =
                bezier_nd::bezier::Bezier::<_, 2, 10>::bernstein_basis_coeff(degree, c, *t);
        }
    }
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

#[track_caller]
fn assert_near_identity(n: usize, m: &[f64]) {
    for i in 0..n * n {
        if i % (n + 1) == 0 {
            assert!((m[i] - 1.0).abs() < 1E-4);
        } else {
            assert!(m[i].abs() < 1E-4);
        }
    }
}

#[track_caller]
fn assert_near_equal(m0: &[f64], m1: &[f64]) {
    for (i, (v0, v1)) in m0.iter().zip(m1.iter()).enumerate() {
        assert!(
            (v0 - v1).abs() < 1E-4,
            "Data {i} is mismatch in {m0:?} <> {m1:?}"
        );
    }
}

#[test]
fn elevate_matrix() {
    let mut m = [0.0f64; 100];

    // Check point to linear matrix
    generate_elevate_by_one_matrix(&mut m, 0);
    assert_eq!(&m[0..2], &[1., 1.]);

    // Check linear - quadratic matrix
    generate_elevate_by_one_matrix(&mut m, 1);
    assert_eq!(&m[0..6], &[1., 0., 0.5, 0.5, 0., 1.]);

    // Check cubic to order 4 matrix
    generate_elevate_by_one_matrix(&mut m, 3);
    assert_eq!(&m[0..4], &[1., 0., 0., 0.]);
    assert_eq!(&m[4..8], &[0.25, 0.75, 0., 0.]);
    assert_eq!(&m[8..12], &[0., 0.5, 0.5, 0.]);
    assert_eq!(&m[12..16], &[0., 0., 0.75, 0.25]);
    assert_eq!(&m[16..20], &[0., 0., 0., 1.]);

    for i in 1..9 {
        generate_elevate_by_one_matrix(&mut m, i);
        eprint!("{i} : [ ");
        let n = (i + 1) * (i + 2);
        for j in 0..n {
            eprint!("{}, ", m[j]);
        }
        eprintln!("]");
        assert_near_equal(bezier_nd::bezier::ELEVATE_BY_ONE_MATRICES_F64[i], &m[0..n]);
    }
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
    generate_elevate_by_one_matrix(&mut elevate, degree);
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
fn bernstein_reduce_matrix_q_to_c() {
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
    generate_elevate_by_one_matrix(&mut elevate, 2);

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
    generate_bs_reduce_matrix(2, &[0., 0.5, 1.0], &mut reduce, &mut elevated_reduce);
    generate_bs_reduce_matrix(
        3,
        &[0., 1.0 / 3.0, 2.0 / 3.0, 1.0],
        &mut reduce,
        &mut elevated_reduce,
    );
    generate_bs_reduce_matrix(
        4,
        &[0., 1.0 / 4.0, 2.0 / 4.0, 3.0 / 4.0, 1.0],
        &mut reduce,
        &mut elevated_reduce,
    );
    generate_bs_reduce_matrix(
        5,
        &[0., 1.0 / 5.0, 2.0 / 5.0, 3.0 / 5.0, 4.0 / 5.0, 1.0],
        &mut reduce,
        &mut elevated_reduce,
    );
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
    assert!(false);
}
