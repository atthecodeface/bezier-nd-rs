//! This library provides functions that return *Metrics* for Bezier curves or the difference between
//! them and other curves or lines.
//!
use crate::{bernstein_fns, utils, BezierEval, BezierMetric, Float, Num};
use geo_nd::vector;

/// Maximum of the squared length of the control points
///
/// This provides a bound on the maximum distance (squared) from the origin of every point on the Bezier curve
pub fn c_sq<F: Num, const D: usize>(bezier: &[[F; D]]) -> F {
    bezier
        .iter()
        .fold(F::ZERO, |acc, pt| utils::max(acc, vector::length_sq(pt)))
}

/// Maximum of the squared length of the control points mapped by the mapping matrix
///
/// This provides a bound on the maximum distance (squared) from the origin of every point on the Bezier curve
pub fn mapped_c_sq<F: Num, const D: usize>(bezier: &[[F; D]], mapping: &[F]) -> F {
    assert_eq!(mapping.len(), bezier.len() * bezier.len());
    mapping.chunks_exact(bezier.len()).fold(F::ZERO, |acc, m| {
        utils::max(acc, vector::length_sq(&vector::sum_scaled(bezier, m)))
    })
}

/// Total of the squared length of the control points
///
/// This is guaranteed to be at least as large as c_sq
pub fn f_sq<F: Num, const D: usize>(bezier: &[[F; D]]) -> F {
    bezier
        .iter()
        .fold(F::ZERO, |acc, pt| acc + vector::length_sq(pt))
}

/// The L2 norm between two beziers given a step dt
///
/// This is an approximation meant for testing/analysis - do not use if
/// performance is required!
///
/// L2 norm is the integral of the distance squared between each point
/// at the same t. This is a simple estimate using summation.
pub fn dl_sq_est<F: Num, const D: usize>(
    bezier: &[[F; D]],
    other: &[[F; D]],
    num_steps: usize,
) -> Option<F> {
    if bezier.len() != other.len() {
        None
    } else {
        Some(
            utils::float_iter(num_steps).fold(F::ZERO, |acc, t| {
                acc + vector::length_sq(&bernstein_fns::values::vector_between_at(bezier, other, t))
            }) / (num_steps as f32).into(),
        )
    }
}

/// The maximum difference between two beziers given a step dt
///
/// This is an approximation meant for testing/analysis - do not use if
/// performance is required!
///
/// Note that dm_sq >= l2, since l2, as can be intuited by it being 1/n.Sum(n)(d_sq(i/n)), for
/// n approaching infinity, and dm_sq is Max(d_sq) = 1/n.Sum(n)(Max(d_sq)), and d_sq() at some t
/// is always at at most Max(d_sq()).
pub fn dm_sq_est<F: Num, const D: usize>(
    bezier: &[[F; D]],
    other: &[[F; D]],
    num_steps: usize,
) -> Option<F> {
    if bezier.len() != other.len() {
        None
    } else {
        Some(utils::float_iter(num_steps).fold(F::ZERO, |acc, t| {
            utils::max(
                acc,
                vector::length_sq(&bernstein_fns::values::vector_between_at(bezier, other, t)),
            )
        }))
    }
}

/// Maximum of the squared length of the difference in control points between two Bezier curves
///
/// This requires the two Bezier curves to be of the same degree, else `None` is returned
///
/// This metric is `Max(P0[j]-P1[j])^2 * (Sum(B[i](t)))^2`, which is
/// the same as `Sum(B[i](t) * Max(P0[i]-P1[i]))^2`; this is more than
/// `Sum(Max(B[i](t) * (P0[i]-P1[i])))^2`
///
/// Hence this value is guaranteed to be greater than or equal to the `dm_sq` metric
///
/// This provides a simply computable bound on the maximum distance (squared) between the respective points on the two Bezier
/// curves for the same parameter `t`.
///
pub fn dc_sq<F: Num, const D: usize>(bezier: &[[F; D]], other: &[[F; D]]) -> Option<F> {
    (bezier.len() == other.len()).then(|| {
        bezier
            .iter()
            .zip(other.iter())
            .fold(F::ZERO, |acc, (b, o)| {
                utils::max(acc, vector::distance_sq(b, o))
            })
    })
}

/// Total of the squared length of the difference in control points between two Bezier curves
///
/// This requires the two Bezier curves to be of the same degree, else `None` is returned
///
/// This is clearly more than the `dc_sq`, the *maximum* of the squared difference between the
/// control points, as that is only one of the `N` control points
///
/// This, though, is clearly less than or equal to `N` times `dc_sq`, the *maximum* of the squared
/// difference between the control points.
pub fn df_sq<F: Num, const D: usize>(bezier: &[[F; D]], other: &[[F; D]]) -> Option<F> {
    (bezier.len() == other.len()).then(|| {
        bezier
            .iter()
            .zip(other.iter())
            .fold(F::ZERO, |acc, (b, o)| acc + vector::distance_sq(b, o))
    })
}

/// Maximum of the squared length of the difference in control points between a bezier and the straight line between its endpoints
///
/// This provides a metric of how far from a straight line the curve is; all points on the
/// Bezier curve will be no further from the straight line than this distance squared
pub fn dc_sq_from_line<F: Num, const D: usize>(bezier: &[[F; D]]) -> F {
    assert!(bezier.len() > 1);
    let l0 = bezier.first().unwrap();
    let l1 = bezier.last().unwrap();
    bezier
        .iter()
        .zip(utils::float_iter(bezier.len()))
        .fold(F::ZERO, |acc, (b, t)| {
            let l = vector::mix(l0, l1, t);
            utils::max(acc, vector::distance_sq(b, &l))
        })
}

/// Total of the squared length of the difference in control points between a bezier and the straight line between its endpoints
///
/// This is guaranteed to be at least as large as `dc_sq_from_line`
pub fn df_sq_from_line<F: Num, const D: usize>(bezier: &[[F; D]]) -> F {
    assert!(bezier.len() > 1);
    let l0 = bezier.first().unwrap();
    let l1 = bezier.last().unwrap();
    bezier
        .iter()
        .zip(utils::float_iter(bezier.len()))
        .fold(F::ZERO, |acc, (b, t)| {
            let l = vector::mix(l0, l1, t);
            acc + vector::distance_sq(b, &l)
        })
}

pub fn metric_from<F: Num, const D: usize>(
    bezier: &[[F; D]],
    other: &[[F; D]],
    metric: BezierMetric,
) -> Option<F> {
    if other.len() != bezier.len() {
        None
    } else if bezier.len() == 0 {
        Some(F::ZERO)
    } else {
        match metric {
            BezierMetric::MaxDistanceSquared(num_steps) => dm_sq_est(bezier, other, num_steps),
            BezierMetric::SumDistanceSquared(num_steps) => dl_sq_est(bezier, other, num_steps),
            BezierMetric::MaxControlSquared => dc_sq(bezier, other),
            _ => df_sq(bezier, other),
        }
    }
}

pub fn metric_from_line<F: Num, const D: usize>(bezier: &[[F; D]], metric: BezierMetric) -> F {
    if bezier.len() <= 2 {
        F::ZERO
    } else {
        match metric {
            BezierMetric::MaxDistanceSquared(_num_steps) => F::ZERO,
            BezierMetric::SumDistanceSquared(_num_steps) => F::ZERO,
            BezierMetric::MaxControlSquared => dc_sq_from_line(bezier),
            _ => df_sq_from_line(bezier),
        }
    }
}

/// Calculates the length of the Bezier when it is rendered down
/// to the given a straightness
///
/// `straightness` is independent of the length of the Bezier
pub fn length_of_lines<F: Float, const D: usize, I: Iterator<Item = ([F; D], [F; D])>>(
    lines: I,
) -> F {
    lines.fold(F::ZERO, |acc, (p0, p1)| {
        acc + geo_nd::vector::distance(&p0, &p1)
    })
}

/// Find the closest point and parameter t on a Bezier to a give point
pub fn t_dsq_closest_to_pt<F: Num, const D: usize>(pts: &[[F; D]], pt: &[F; D]) -> Option<(F, F)> {
    match pts.len() {
        1 => Some((F::ZERO, vector::distance_sq(&pts[0], pt))),
        2 => <&[[F; D]] as TryInto<&[[F; D]; 2]>>::try_into(&pts[0..2])
            .unwrap()
            .t_dsq_closest_to_pt(pt),
        3 => <&[[F; D]] as TryInto<&[[F; D]; 3]>>::try_into(&pts[0..3])
            .unwrap()
            .t_dsq_closest_to_pt(pt),
        4 => <&[[F; D]] as TryInto<&[[F; D]; 4]>>::try_into(&pts[0..4])
            .unwrap()
            .t_dsq_closest_to_pt(pt),
        _ => None,
    }
}

/// Estimate the minimum squared distance to a Bezier
pub fn est_min_distance_sq_to<F: Num, const D: usize>(pts: &[[F; D]], pt: &[F; D]) -> F {
    // eprintln!("est_min_distance_sq_to {pts:?} {pt:?}");
    let l0 = pts.first().unwrap();
    let l1 = pts.last().unwrap();
    let d_sq = crate::utils::distance_sq_to_line_segment(pt, l0, l1);
    let t_iter = utils::float_iter(pts.len());
    let dc_sq = pts.iter().zip(t_iter).fold(F::ZERO, |acc, (pt, t)| {
        utils::max(acc, vector::distance_sq(pt, &vector::mix(l0, l1, t)))
    });
    // min possible distance to Bezier = sqrt(d_sq) - sqrt(dc_sq)
    if d_sq < dc_sq {
        F::ZERO
    } else {
        // dout <= min possible distance = sqrt(d_sq) - sqrt(dc_sq)
        utils::est_d_m_c_from_dsq_m_dcsq(d_sq, dc_sq)
    }
}
