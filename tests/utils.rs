use bezier_nd::Bezier;
use geo_nd::Float;
use geo_nd::{
    matrix,
    vector::{self, reduce},
    FArray, Num,
};
use num::rational::Rational64;

pub fn float_iter<N: Num>(t0: N, t1: N, n: usize) -> impl Iterator<Item = N> {
    assert!(
        n >= 2,
        "Float iterator must have at least two steps for begin and end"
    );
    let r = t1 - t0;
    assert!(r >= N::zero(), "Float range must be positive");
    (0..n).map(move |i: usize| (r * N::from_usize(i).unwrap() / N::from_usize(n - 1).unwrap()))
}

#[track_caller]
pub fn assert_near_identity<N: Num>(n: usize, m: &[N]) {
    let eps = N::from_usize(1).unwrap() / N::from_usize(10000).unwrap();
    for i in 0..n * n {
        let exp = {
            if i % (n + 1) == 0 {
                N::one()
            } else {
                N::zero()
            }
        };
        assert!(m[i] - exp > -eps && m[i] - exp < eps);
    }
}

#[track_caller]
pub fn assert_near_equal<N: Num>(m0: &[N], m1: &[N]) {
    let eps = N::from_usize(1).unwrap() / N::from_usize(10000).unwrap();

    for (i, (v0, v1)) in m0.iter().zip(m1.iter()).enumerate() {
        let dv = *v0 - *v1;
        assert!(
            dv >= -eps && dv <= eps,
            "Data {i} is mismatch in {m0:?} <> {m1:?}"
        );
    }
}

#[track_caller]
pub fn assert_near_equal_scale<N: Num>(m0: &[N], m1: &[N], scale: N) {
    let eps = N::from_usize(1).unwrap() / N::from_usize(10000).unwrap();

    for (i, (v0, v1)) in m0.iter().zip(m1.iter()).enumerate() {
        let dv = (*v0) * scale - *v1;
        assert!(
            dv >= -eps && dv <= eps,
            "Data {i} is mismatch in {m0:?} <> {m1:?}"
        );
    }
}

pub fn test_subsection(b: &Bezier<f64, 2>, sub: &Bezier<f64, 2>, t0: f64, t1: f64) {
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

pub fn test_beziers_approx_eq(b0: &Bezier<f64, 2>, b1: &Bezier<f64, 2>) {
    test_subsection(b0, b1, 0.0, 1.0);
    test_subsection(b0, &b1.bezier_between(0.0, 1.0), 0.0, 1.0);
    test_subsection(b0, &b1.bezier_between(0.1, 0.4), 0.1, 0.4);
    test_subsection(b0, &b1.bisect().0, 0.0, 0.5);
    test_subsection(b0, &b1.bisect().1, 0.5, 1.0);
}

#[track_caller]
pub fn generate_bernstein_matrix<F: Float>(matrix: &mut [F], degree: usize, ts: &[F]) {
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
                bezier_nd::bernstein::bezier_fns::basis_coeff(degree, c, *t);
        }
    }
}

pub fn bernstein_basis_coeff_br(degree: usize, i: usize, t: Rational64) -> Rational64 {
    let u = Rational64::ONE - t;
    let mut result = Rational64::ONE;
    for c in 0..degree {
        if c < i {
            result *= t;
        } else {
            result *= u;
        }
    }
    // Multiply by n! / (n-c)!
    for j in 0..i {
        eprintln!("Mult by {} {degree} {j} {i}", degree - j);
        let f: Rational64 = ((degree - j) as i64).into();
        result *= f;
    }
    // Divide by c!
    for j in 1..=i {
        eprintln!("Divide by {}", j);
        let f: Rational64 = (j as i64).into();
        result /= f;
    }

    dbg!(&degree, &i, &t, &u, &result);
    result
}

#[track_caller]
pub fn generate_bernstein_matrix_br(matrix: &mut [Rational64], degree: usize, ts: &[Rational64]) {
    assert_eq!(
        matrix.len(),
        (degree + 1) * ts.len(),
        "Trying to generate an {} by {} matrix, so must be given a matrix of that size",
        ts.len(),
        degree + 1
    );
    for (r, t) in ts.iter().enumerate() {
        for c in 0..(degree + 1) {
            matrix[r * (degree + 1) + c] = bernstein_basis_coeff_br(degree, c, *t);
        }
    }
}
