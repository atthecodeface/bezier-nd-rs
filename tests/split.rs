//a Imports
use bezier_nd::Bezier;
use geo_nd::{matrix, vector, FArray};

fn float_iter(t0: f32, t1: f32, n: usize) -> impl Iterator<Item = f32> {
    assert!(
        n >= 2,
        "Float iterator must have at least two steps for begin and end"
    );
    let r = t1 - t0;
    assert!(r >= 0.0, "Float range must be positive");
    (0..n).map(move |i: usize| (i as f32) / (r * ((n - 1) as f32)))
}

fn test_subsection(b: &Bezier<f32, 2>, sub: &Bezier<f32, 2>, t0: f32, t1: f32) {
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

fn test_beziers_approx_eq(b0: &Bezier<f32, 2>, b1: &Bezier<f32, 2>) {
    test_subsection(b0, b1, 0.0, 1.0);
    test_subsection(b0, &b1.bezier_between(0.0, 1.0), 0.0, 1.0);
    test_subsection(b0, &b1.bezier_between(0.1, 0.4), 0.1, 0.4);
    test_subsection(b0, &b1.bisect().0, 0.0, 0.5);
    test_subsection(b0, &b1.bisect().1, 0.5, 1.0);
}

#[test]
fn bisect() {
    let p0: FArray<f32, 2> = [0., 0.].into();
    let p1: FArray<f32, 2> = [10., 0.].into();
    let p2: FArray<f32, 2> = [6., 1.].into();
    let p3: FArray<f32, 2> = [20., 5.].into();

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
    let p0: FArray<f32, 2> = [0., 0.].into();
    let p1: FArray<f32, 2> = [10., 0.].into();
    let p2: FArray<f32, 2> = [6., 1.].into();
    let p3: FArray<f32, 2> = [20., 5.].into();

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

fn generate_elevate_by_one_matrix(matrix: &mut [f32], degree: usize) {
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

fn generate_bernstein_matrix(matrix: &mut [f32], degree: usize, ts: &[f32]) {
    assert_eq!(
        matrix.len(),
        (degree + 1) * ts.len(),
        "Trying to generate an {} by {} matrix, so must be given a matrix of that size",
        ts.len(),
        degree + 1
    );
    for (r, t) in ts.iter().enumerate() {
        for c in 0..(degree + 1) {
            // matrix[r * (degree + 1) + c] =
            matrix[c * ts.len() + r] =
                bezier_nd::bezier::Bezier::<_, 2, 10>::bernstein_basis_coeff(degree, c, *t);
        }
    }
}

#[test]
fn bernstein_matrix() {
    let mut bern_n = [0.0f32; 100];
    generate_bernstein_matrix(&mut bern_n[0..3], 2, &[0.]);
    assert_eq!(&bern_n[0..3], &[1., 0., 0.]);
    generate_bernstein_matrix(&mut bern_n[0..3], 2, &[0.5]);
    assert_eq!(&bern_n[0..3], &[0.25, 0.5, 0.25]);
    generate_bernstein_matrix(&mut bern_n[0..3], 2, &[1.0]);
    assert_eq!(&bern_n[0..3], &[0., 0., 1.]);
}

#[track_caller]
fn assert_near_equal(m0: &[f32], m1: &[f32]) {
    for (i, (v0, v1)) in m0.iter().zip(m1.iter()).enumerate() {
        assert!(
            (v0 - v1).abs() < 1E-4,
            "Data {i} is mismatch in {m0:?} <> {m1:?}"
        );
    }
}

#[test]
fn bernstein_reduce_matrix() {
    let mut bern_n = [0.0f32; 100];
    let mut bern_np1 = [0.0f32; 100];
    let ts = [0., 0.5, 1.0];
    generate_bernstein_matrix(&mut bern_n[0..9], 2, &ts);
    assert_near_equal(&bern_n, &[1.0, 0.25, 0.0, 0., 0.5, 0., 0., 0.25, 1.]);
    generate_bernstein_matrix(&mut bern_np1[0..12], 3, &ts);
    assert_near_equal(
        &bern_np1,
        &[1.0, 0.125, 0., 0., 0.375, 0., 0., 0.375, 0., 0., 0.125, 1.],
    );
    let mut bern_n_inverse = bern_n.clone();
    assert!(
        matrix::invert::<_, 9, 3>(&mut bern_n_inverse.first_chunk_mut::<9>().unwrap()),
        "Bernstein basis must be invertible"
    );
    assert_near_equal(
        &bern_n_inverse,
        &[1.0, -0.5, 0.0, 0., 2.0, 0., 0., -0.5, 1.],
    );
    let reduce = matrix::multiply::<_, 12, 9, 12, 4, 3, 3>(
        bern_np1.first_chunk::<12>().unwrap(),
        bern_n_inverse.first_chunk::<9>().unwrap(),
    );
    eprintln!("bernstein_reduce[3] = {reduce:?}");
    assert_near_equal(
        &reduce,
        &[
            1.0, -0.25, 0.0, 0.0, 0.75, 0.0, 0.0, 0.75, 0.0, 0.0, -0.25, 1.0,
        ],
    );
}

#[test]
fn elevate_matrix() {
    let mut m = [0.0f32; 100];

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
        assert_eq!(bezier_nd::bezier::ELEVATE_BY_ONE_MATRICES[i], &m[0..n]);
    }
}
