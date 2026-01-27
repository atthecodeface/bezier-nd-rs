use bezier_nd::BezierEval;
use bezier_nd::Num;
use geo_nd::vector;

use super::float_iter;

/// Find the maximum distance between the Bezier for 0<=t<=1 in 'steps' intervals
pub fn max_distance_sq<
    F: Num,
    const D: usize,
    B0: BezierEval<F, [F; D]>,
    B1: BezierEval<F, [F; D]>,
>(
    b0: &B0,
    b1: &B1,
    steps: usize,
) -> F {
    let mut max_d_sq = F::ZERO;
    for t in float_iter(F::ZERO, F::ONE, steps) {
        let p0 = b0.point_at(t);
        let p1 = b1.point_at(t);
        let d_sq = vector::distance_sq(&p0, &p1);
        eprintln!("t: {t} pts: {p0:?} {p1:?} d_sq:{d_sq}");
        if d_sq > max_d_sq {
            max_d_sq = d_sq;
        }
    }
    max_d_sq
}
