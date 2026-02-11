use crate::Num;
use geo_nd::vector;

/// Calculate a point on Bernstein Bezier at parameter t
pub fn point_at<F: Num, const D: usize>(pts: &[[F; D]], t: F) -> [F; D] {
    let degree = pts.len() - 1;
    pts.iter()
        .zip(super::basis_coeff_enum_num(degree, t))
        .fold([F::ZERO; D], |acc, (pt, (_i, coeff))| {
            vector::add(acc, pt, coeff)
        })
}

/// Calculate the first derivative of a Bernstein Bezier at parameter t
pub fn derivative_at<F: Num, const D: usize>(pts: &[[F; D]], t: F) -> (F, [F; D]) {
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
