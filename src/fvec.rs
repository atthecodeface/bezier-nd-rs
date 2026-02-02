use crate::bernstein::bezier_fns;
use crate::utils;
use crate::{
    BezierDistance, BezierEval, BezierMinMax, BezierReduce, BezierSection, BezierSplit, BoxedBezier,
};
use crate::{Float, Num};
use geo_nd::vector;

impl<F: Float, const D: usize> BezierEval<F, [F; D]> for Vec<[F; D]> {
    fn point_at(&self, t: F) -> [F; D] {
        let degree = self.len() - 1;

        self.iter()
            .zip(bezier_fns::basis_coeff_enum(degree, t))
            .fold([F::ZERO; D], |acc, (pt, (_i, coeff))| {
                vector::add(acc, pt, coeff)
            })
    }
    fn derivative_at(&self, _t: F) -> (F, [F; D]) {
        (F::ONE, vector::add(self[1], &self[0], -F::ONE))
    }
    fn endpoints(&self) -> (&[F; D], &[F; D]) {
        (&self[0], &self[1])
    }
    fn closeness_sq_to_line(&self) -> F {
        F::ZERO
    }
    fn dc_sq_from_line(&self) -> F {
        F::ZERO
    }

    fn num_control_points(&self) -> usize {
        2
    }
    fn control_point(&self, n: usize) -> &[F; D] {
        &self[n]
    }
}
