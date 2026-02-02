use crate::utils;
use crate::{
    bernstein_fns, BezierDistance, BezierEval, BezierMinMax, BezierReduce, BezierSection,
    BezierSplit, BoxedBezier,
};
use crate::{Float, Num};
use geo_nd::vector;

impl<F: Float, const D: usize> BezierEval<F, [F; D]> for Vec<[F; D]> {
    fn point_at(&self, t: F) -> [F; D] {
        let degree = self.len() - 1;

        self.iter()
            .zip(bernstein_fns::basis_coeff_enum(degree, t))
            .fold([F::ZERO; D], |acc, (pt, (_i, coeff))| {
                vector::add(acc, pt, coeff)
            })
    }
    fn derivative_at(&self, t: F) -> (F, [F; D]) {
        let degree = self.len() - 1;
        let (reduce, coeffs) = bernstein_fns::basis_dt_coeff_enum(degree, t);
        (
            reduce,
            self.iter()
                .zip(coeffs)
                .fold([F::ZERO; D], |acc, (pt, (_i, coeff))| {
                    vector::add(acc, pt, coeff)
                }),
        )
    }
    fn endpoints(&self) -> (&[F; D], &[F; D]) {
        (self.first().unwrap(), self.last().unwrap())
    }
    fn closeness_sq_to_line(&self) -> F {
        self.dc_sq_from_line()
    }
    fn dc_sq_from_line(&self) -> F {
        bernstein_fns::metrics::dc_sq_from_line(
            self.iter(),
            self.first().unwrap(),
            self.last().unwrap(),
        )
    }

    fn num_control_points(&self) -> usize {
        self.len()
    }
    fn control_point(&self, n: usize) -> &[F; D] {
        &self[n]
    }
}

impl<F: Num, const D: usize> BezierSplit for Vec<[F; D]> {
    fn split(&self) -> (Self, Self) {
        let mut to_split = self.clone();
        let mut b0 = self.clone();
        let mut b1 = self.clone();
        bernstein_fns::split::into_two_at_de_cast(&mut to_split, 0.5_f32.into(), &mut b0, &mut b1);
        (b0, b1)
    }
}

impl<F: Num, const D: usize> BezierSection<F> for Vec<[F; D]> {
    /// Split the Bezier into two (one covering parameter values 0..t, the other t..1)
    fn split_at(&self, t: F) -> (Self, Self) {
        let mut to_split = self.clone();
        let mut b0 = self.clone();
        let mut b1 = self.clone();
        bernstein_fns::split::into_two_at_de_cast(&mut to_split, t, &mut b0, &mut b1);
        (b0, b1)
    }

    /// Generate the Bezier with parameter t' where 0<=t'<=1 provides the same points
    /// as the original Bezier with t0<=t<=1
    fn section(&self, t0: F, t1: F) -> Self {
        let mut to_split = self.clone();
        let mut b0 = self.clone();
        let mut b1 = self.clone();
        let t10 = (t1 - t0) / (F::ONE - t0);
        bernstein_fns::split::into_two_at_de_cast(&mut to_split, t0, &mut b0, &mut b1);
        bernstein_fns::split::into_two_at_de_cast(&mut b1, t10, &mut to_split, &mut b0);
        to_split
    }
}
