use crate::{utils, Num};
use geo_nd::vector;

/// Maximum of the squared length of the control points
pub fn c_sq<'a, F: Num, const D: usize, I: Iterator<Item = &'a [F; D]>>(pts: I) -> F {
    pts.fold(F::ZERO, |acc, pt| utils::max(acc, vector::length_sq(pt)))
}

/// Total of the squared length of the control points
pub fn f_sq<'a, F: Num, const D: usize, I: Iterator<Item = &'a [F; D]>>(pts: I) -> F {
    pts.fold(F::ZERO, |acc, pt| acc + vector::length_sq(pt))
}

/// Maximum of the squared length of the control points
pub fn dc_sq_from_line<'a, F: Num, const D: usize, I: ExactSizeIterator<Item = &'a [F; D]> + 'a>(
    pts: I,
    l0: &[F; D],
    l1: &[F; D],
) -> F {
    let n = pts.len();
    pts.zip(utils::float_iter(F::ZERO, F::ONE, n))
        .fold(F::ZERO, |acc, (pt, t)| {
            utils::max(acc, vector::distance_sq(pt, &vector::mix(l0, l1, t)))
        })
}

/// Total of the squared length of the control points
pub fn df_sq_from_line<'a, F: Num, const D: usize, I: ExactSizeIterator<Item = &'a [F; D]> + 'a>(
    pts: I,
    l0: &[F; D],
    l1: &[F; D],
) -> F {
    let n = pts.len();
    pts.zip(utils::float_iter(F::ZERO, F::ONE, n))
        .fold(F::ZERO, |acc, (pt, t)| {
            acc + vector::distance_sq(pt, &vector::mix(l0, l1, t))
        })
}
