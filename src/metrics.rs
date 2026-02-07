use crate::{utils, BezierEval, Num};
use geo_nd::vector;

pub fn c_sq<F: Num, const D: usize, B: BezierEval<F, [F; D]>>(bezier: &B) -> F {
    let mut result = F::ZERO;
    bezier.for_each_control_point(&mut |_i, pt| {
        result = utils::max(result, vector::length_sq(pt));
    });
    result
}

pub fn f_sq<F: Num, const D: usize, B: BezierEval<F, [F; D]>>(bezier: &B) -> F {
    let mut result = F::ZERO;
    bezier.for_each_control_point(&mut |_i, pt| {
        result += vector::length_sq(pt);
    });
    result
}

pub fn dc_sq<F: Num, const D: usize, B: BezierEval<F, [F; D]>>(bezier: &B, other: &B) -> Option<F> {
    if bezier.num_control_points() != other.num_control_points() {
        None
    } else {
        let mut result = F::ZERO;
        bezier.for_each_control_point(&mut |i, pt| {
            result = utils::max(result, vector::distance_sq(pt, other.control_point(i)));
        });
        Some(result)
    }
}

pub fn df_sq<F: Num, const D: usize, B: BezierEval<F, [F; D]>>(bezier: &B, other: &B) -> Option<F> {
    if bezier.num_control_points() != other.num_control_points() {
        None
    } else {
        let mut result = F::ZERO;
        bezier.for_each_control_point(&mut |i, pt| {
            result += vector::distance_sq(pt, other.control_point(i));
        });
        Some(result)
    }
}

pub fn dc_sq_from_line<F: Num, const D: usize, B: BezierEval<F, [F; D]>>(bezier: &B) -> F {
    let mut result = F::ZERO;
    let mut iter = utils::float_iter(bezier.num_control_points());
    let (l0, l1) = bezier.endpoints();
    bezier.for_each_control_point(&mut |_i, pt| {
        let l = vector::mix(l0, l1, iter.next().unwrap());
        result = utils::max(result, vector::distance_sq(pt, &l));
    });
    result
}

pub fn df_sq_from_line<F: Num, const D: usize, B: BezierEval<F, [F; D]>>(bezier: &B) -> F {
    let mut result = F::ZERO;
    let mut iter = utils::float_iter(bezier.num_control_points());
    let (l0, l1) = bezier.endpoints();
    bezier.for_each_control_point(&mut |_i, pt| {
        let l = vector::mix(l0, l1, iter.next().unwrap());
        result += vector::distance_sq(pt, &l);
    });
    result
}
