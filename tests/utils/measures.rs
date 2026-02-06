use bezier_nd::BezierEval;
use bezier_nd::Num;
use geo_nd::vector;

use super::{float_iter, float_iter_between};

/// Find the maximum distance between the Bezier for 0<=t<=1 in 'steps' intervals
#[must_use]
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
    for t in float_iter(steps) {
        let p0 = b0.point_at(t);
        let p1 = b1.point_at(t);
        let d_sq = vector::distance_sq(&p0, &p1);
        // eprintln!("t: {t} pts: {p0:?} {p1:?} d_sq:{d_sq}");
        if d_sq > max_d_sq {
            max_d_sq = d_sq;
        }
    }
    max_d_sq
}

/// Find the maximum distance between `b0` for t0<=t<=t1 and `b1` for 0<=t<=1 in 'steps' intervals
#[must_use]
pub fn max_distance_sq_subsection<
    F: Num,
    const D: usize,
    B0: BezierEval<F, [F; D]>,
    B1: BezierEval<F, [F; D]>,
>(
    b0: &B0,
    b1: &B1,
    t0: F,
    t1: F,
    steps: usize,
) -> F {
    let mut max_d_sq = F::ZERO;
    for (b0_t, b1_t) in float_iter_between(t0, t1, steps).zip(float_iter(steps)) {
        let p0 = b0.point_at(b0_t);
        let p1 = b1.point_at(b1_t);
        let d_sq = vector::distance_sq(&p0, &p1);
        // eprintln!("t: {b0_t} <> {b1_t} pts: {p0:?} {p1:?} d_sq:{d_sq}");
        if d_sq > max_d_sq {
            max_d_sq = d_sq;
        }
    }
    max_d_sq
}

pub fn assert_pts_all_within_closeness_sq_of_pts<
    'a,
    F: Num,
    I: Iterator<Item = &'a [F; D]> + 'a,
    const D: usize,
>(
    good_pts: &[[F; D]],
    test_pts: I,
    closeness_sq: F,
) {
    let mut s = 0;
    let ngood = good_pts.len();
    for test_pt in test_pts {
        let mut found = false;
        for t in 0..(2 * ngood) {
            let i = {
                if (t & 1) != 0 {
                    if s < 1 + t / 2 {
                        continue;
                    } else {
                        s - 1 + t / 2
                    }
                } else {
                    s + t / 2
                }
            };
            if i >= ngood {
                continue;
            }
            let d_sq = vector::distance_sq(&test_pt, &good_pts[i]);
            if d_sq <= closeness_sq {
                s = i;
                found = true;
                break;
            }
        }
        assert!(
            found,
            "Failed to find a point in segment_pts within {closeness_sq} of test point {test_pt:?}"
        );
    }
}
