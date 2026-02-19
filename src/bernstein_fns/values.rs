use crate::Num;
use geo_nd::vector;

/// Calculate a point on Bernstein Bezier at parameter t
#[track_caller]
pub fn point_at<F: Num, const D: usize>(pts: &[[F; D]], t: F) -> [F; D] {
    assert!(!pts.is_empty());
    let degree = pts.len() - 1;
    pts.iter()
        .zip(super::basis_coeff_enum_num(degree, t))
        .fold([F::ZERO; D], |acc, (pt, (_i, coeff))| {
            vector::add(acc, pt, coeff)
        })
}

/// Calculate the vector between points on two Bezier curves with the same `t`
#[track_caller]
pub fn vector_between_at<F: Num, const D: usize>(pts: &[[F; D]], other: &[[F; D]], t: F) -> [F; D] {
    assert!(!pts.is_empty());
    let degree = pts.len() - 1;
    super::basis_coeff_enum_num(degree, t)
        .zip(pts.iter().zip(other.iter()))
        .fold([F::ZERO; D], |acc, ((_i, coeff), (p, o))| {
            vector::sum_scaled(&[acc, *p, *o], &[F::ONE, coeff, -coeff])
        })
}

/// Calculate the first derivative of a Bernstein Bezier at parameter t
#[track_caller]
pub fn derivative_at<F: Num, const D: usize>(pts: &[[F; D]], t: F) -> (F, [F; D]) {
    assert!(!pts.is_empty());
    let degree = pts.len() - 1;
    let (reduce, coeffs) = super::basis_dt_coeff_enum_num(degree, t);
    (
        reduce,
        pts.iter()
            .zip(coeffs)
            .fold([F::ZERO; D], |acc, (pt, (_i, coeff))| {
                vector::add(acc, pt, coeff)
            }),
    )
}

/// Determine which side of a line a point is in 2D
pub fn side_of_line<F: Num>(pt: &[F; 2], p0: &[F; 2], p1: &[F; 2]) -> F {
    (pt[0] - p0[0]) * (p1[1] - p0[1]) - (pt[1] - p0[1]) * (p1[0] - p0[0])
}
