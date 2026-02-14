use crate::Num;
use geo_nd::vector;

/// Transform a slice of points by a matrix into a (new) slice of points
pub fn transform_pts<F: Num, const D: usize>(matrix: &[F], pts: &[[F; D]], result: &mut [[F; D]]) {
    assert_eq!(
        pts.len() * result.len(),
        matrix.len(),
        "Transformation matrix must map {} pts to {} pts",
        pts.len(),
        result.len()
    );
    for (r, m) in result.iter_mut().zip(matrix.chunks_exact(pts.len())) {
        *r = vector::sum_scaled(pts, m);
    }
}
