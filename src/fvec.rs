use crate::{bernstein_fns, BezierEval, BezierSection, BezierSplit};
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
        self.split_at(0.5_f32.into())
    }
}

impl<F: Num, const D: usize> BezierSection<F> for Vec<[F; D]> {
    /// Split the Bezier into two (one covering parameter values 0..t, the other t..1)
    fn split_at(&self, t: F) -> (Self, Self) {
        let mut latter = self.clone();
        let mut first = self.clone();
        bernstein_fns::split::into_two_at_de_cast(&mut latter, t, &mut first);
        (first, latter)
    }

    /// Generate the Bezier with parameter t' where 0<=t'<=1 provides the same points
    /// as the original Bezier with t0<=t<=1
    fn section(&self, t0: F, t1: F) -> Self {
        let mut to_split = self.clone();
        if t0 > F::ZERO {
            bernstein_fns::split::bezier_from_de_cast(&mut to_split, t0);
        }
        if t1 < F::ONE {
            let t10 = (t1 - t0) / (F::ONE - t0);
            bernstein_fns::split::bezier_to_de_cast(&mut to_split, t10);
        }
        to_split
    }
}
