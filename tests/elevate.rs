mod utils;
use bezier_nd::Bezier;
use geo_nd::FArray;
use utils::test_beziers_approx_eq;
use utils::assert_near_equal_scale;

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

#[test]
fn elevate_matrix_isize() {
    let mut m = [0_isize; 100];

    // Check point to linear matrix
    assert_eq!(
        bezier_nd::bernstein::bezier_fns::generate_elevate_by_one_matrix(&mut m, 0),
        1
    );
    assert_eq!(&m[0..2], &[1, 1]);

    // Check linear - quadratic matrix
    assert_eq!(
        bezier_nd::bernstein::bezier_fns::generate_elevate_by_one_matrix(&mut m, 1),
        2
    );
    assert_eq!(&m[0..6], &[2, 0, 1, 1, 0, 2]);

    // Check cubic to order 4 matrix
    assert_eq!(
        bezier_nd::bernstein::bezier_fns::generate_elevate_by_one_matrix(&mut m, 3),
        4
    );
    assert_eq!(&m[0..4], &[4, 0, 0, 0]);
    assert_eq!(&m[4..8], &[1, 3, 0, 0]);
    assert_eq!(&m[8..12], &[0, 2, 2, 0]);
    assert_eq!(&m[12..16], &[0, 0, 3, 1]);
    assert_eq!(&m[16..20], &[0, 0, 0, 4]);

    for i in 1..9 {
        let scale = bezier_nd::bernstein::bezier_fns::generate_elevate_by_one_matrix(&mut m, i);
        eprint!("&[ ");
        let n = (i + 1) * (i + 2);
        for j in 0..n {
            eprint!("{}., ", m[j]);
        }
        eprintln!("],");
        let m_f64: Vec<f64> = m.iter().map(|m| (*m) as f64).collect();
        assert_near_equal_scale(bezier_nd::ELEVATE_BY_ONE_MATRICES_F64[i], &m_f64[0..n], 1.0);
    }
}

#[test]
fn elevate_matrix() {
    let mut m = [0.0f64; 100];

    // Check point to linear matrix
    assert_eq!(
        bezier_nd::bernstein::bezier_fns::generate_elevate_by_one_matrix(&mut m, 0),
        1.0
    );
    assert_eq!(&m[0..2], &[1., 1.]);

    // Check linear - quadratic matrix
    assert_eq!(
        bezier_nd::bernstein::bezier_fns::generate_elevate_by_one_matrix(&mut m, 1),
        2.0
    );
    assert_eq!(&m[0..6], &[2., 0., 1.0, 1.0, 0., 2.]);

    // Check cubic to order 4 matrix
    assert_eq!(
        bezier_nd::bernstein::bezier_fns::generate_elevate_by_one_matrix(&mut m, 3),
        4.0
    );
    assert_eq!(&m[0..4], &[4., 0., 0., 0.]);
    assert_eq!(&m[4..8], &[1.0, 3., 0., 0.]);
    assert_eq!(&m[8..12], &[0., 2., 2., 0.]);
    assert_eq!(&m[12..16], &[0., 0., 3., 1.]);
    assert_eq!(&m[16..20], &[0., 0., 0., 4.]);
}
