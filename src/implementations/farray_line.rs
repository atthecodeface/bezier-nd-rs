use crate::Num;
use crate::{
    metrics, utils, BezierBuilder, BezierConstruct, BezierElevate, BezierError, BezierEval,
    BezierFlatIterator, BezierIterationType, BezierMetric, BezierOps, BezierReduce,
    BezierReduction, BezierSplit, BoxedBezier,
};

use geo_nd::vector;

impl<F: Num, const D: usize> BezierEval<F, [F; D]> for [[F; D]; 2] {
    fn distance_sq_between(&self, p0: &[F; D], p1: &[F; D]) -> F {
        vector::distance_sq(p0, p1)
    }
    fn point_at(&self, t: F) -> [F; D] {
        vector::add(vector::scale(self[0], F::ONE - t), &self[1], t)
    }
    fn derivative_at(&self, _t: F) -> (F, [F; D]) {
        (F::ONE, vector::add(self[1], &self[0], -F::ONE))
    }
    fn closeness_sq_to_line(&self) -> F {
        F::ZERO
    }
    fn dc_sq_from_line(&self) -> F {
        F::ZERO
    }
    fn endpoints(&self) -> ([F; D], [F; D]) {
        (self[0], self[1])
    }
    fn control_points(&self) -> &[[F; D]] {
        self
    }
    fn metric_from(&self, other: Option<&[[F; D]]>, metric: BezierMetric) -> Option<F> {
        if let Some(other) = other {
            metrics::metric_from(self, other, metric)
        } else {
            Some(metrics::metric_from_line(self, metric))
        }
    }

    /// A value of t and distance squared to it, 0<=t<=1, for which the distance between the point
    ///  P and the line between the two points is a minimum.
    ///
    /// The distance between a line (p0, p1) and a point P has a single minimum as one moves along the line.
    ///
    /// For every point Q on the line q, (q-p0) = t(p1-p0) for some t; (q-p0).(p1-p0) = t.|p1-p0|^2,
    /// or t = (q-p0).(p1-p0) / |p1-p0|^2
    ///
    /// If this value returns a value 0<=t<=1 then no other point on the Bezier with t in that range
    /// has a smaller distance to the point P
    ///
    /// If this value returns a value t < 0 then the Bezier at t=0 will be the closest point on the
    /// Bezier with parameter 0<=t<=1 to point P
    ///
    /// If this value returns a value t > 1 then the Bezier at t=1 will be the closest point on the
    /// Bezier with parameter 0<=t<=1 to point P
    ///
    fn t_dsq_closest_to_pt(&self, pt: &[F; D]) -> Option<(F, F)> {
        let (t_times_len_sq, line_len_sq, valid) = utils::relative_to_line(pt, &self[0], &self[1]);
        if valid {
            if t_times_len_sq < F::ZERO {
                Some((F::ZERO, vector::distance_sq(pt, &self[0])))
            } else if t_times_len_sq > line_len_sq {
                Some((F::ONE, vector::distance_sq(pt, &self[1])))
            } else {
                Some((
                    t_times_len_sq / line_len_sq,
                    vector::length_sq(&vector::sub(self[0], pt, F::ONE)) - t_times_len_sq,
                ))
            }
        } else {
            Some((F::ZERO, vector::distance_sq(pt, &self[0])))
        }
    }

    /// An estimate of the minimum distance squared from the Bezier to the point
    /// for 0<=t<=1.
    ///
    /// This returns the *actual* minimum distance to the line segmenr from the point
    fn est_min_distance_sq_to(&self, p: &[F; D]) -> F {
        utils::distance_sq_to_line_segment(p, &self[0], &self[1])
    }

    fn t_coords_at_min_max(
        &self,
        pt_index: usize,
        give_min: bool,
        give_max: bool,
    ) -> (Option<(F, F)>, Option<(F, F)>) {
        utils::opt_min_and_max_tc(
            give_min,
            give_max,
            (F::ZERO, self[0][pt_index]),
            (F::ONE, self[1][pt_index]),
            None,
        )
    }
}

impl<F: Num, const D: usize> BezierOps<F, [F; D]> for [[F; D]; 2] {
    fn add(&mut self, other: &Self) -> bool {
        for (s, o) in self.iter_mut().zip(other.iter()) {
            *s = vector::add(*s, o, F::ONE)
        }
        true
    }
    fn sub(&mut self, other: &Self) -> bool {
        for (s, o) in self.iter_mut().zip(other.iter()) {
            *s = vector::sub(*s, o, F::ONE)
        }
        true
    }
    fn scale(&mut self, scale: F) {
        for s in self.iter_mut() {
            *s = vector::scale(*s, scale);
        }
    }
    fn map_pts(&mut self, map: &dyn Fn(usize, &[F; D]) -> [F; D]) {
        for (i, s) in self.iter_mut().enumerate() {
            *s = map(i, s);
        }
    }
    fn map_all_pts<'a>(&'a mut self, map: &'a mut dyn FnMut(&'a mut [[F; D]]) -> bool) -> bool {
        map(self)
    }
}

impl<F: Num, const D: usize> BezierSplit<F> for [[F; D]; 2] {
    fn split(&self) -> (Self, Self) {
        let m = vector::sum_scaled(self, &[F::frac(1, 2), F::frac(1, 2)]);
        ([self[0], m], [m, self[1]])
    }

    /// Split the Bezier into two (one covering parameter values 0..t, the other t..1)
    fn split_at(&self, t: F) -> (Self, Self) {
        let m = vector::sum_scaled(self, &[F::ONE - t, t]);
        ([self[0], m], [m, self[1]])
    }

    /// Generate the Bezier with parameter t' where 0<=t'<=1 provides the same points
    /// as the original Bezier with t0<=t<=1
    fn section(&self, t0: F, t1: F) -> Self {
        let m0 = vector::sum_scaled(self, &[F::ONE - t0, t0]);
        let m1 = vector::sum_scaled(self, &[F::ONE - t1, t1]);
        [m0, m1]
    }
}
impl<F, const D: usize> BezierFlatIterator<F, [F; D]> for [[F; D]; 2]
where
    F: Num,
{
    /// Return an iterator of N points along the Bezier at even steps of `t`
    fn as_t_lines(
        &self,
        iter_type: BezierIterationType<F>,
    ) -> impl Iterator<Item = (F, [F; D], F, [F; D])> {
        let n = if let BezierIterationType::Uniform(n) = iter_type {
            n
        } else {
            2
        };
        utils::float_iter(n)
            .zip(utils::float_iter(n).skip(1))
            .map(|(t0, t1)| (t0, self.point_at(t0), t1, self.point_at(t1)))
    }

    /// Return an iterator of N points along the Bezier at even steps of `t`
    fn as_t_points(&self, iter_type: BezierIterationType<F>) -> impl Iterator<Item = (F, [F; D])> {
        let n = if let BezierIterationType::Uniform(n) = iter_type {
            n
        } else {
            2
        };
        utils::float_iter(n).map(|t| (t, self.point_at(t)))
    }
}

impl<F: Num, const D: usize> BoxedBezier<F, [F; D]> for [[F; D]; 2] {
    fn boxed_reduce(&self) -> Option<Box<dyn BoxedBezier<F, [F; D]>>> {
        None
    }
    fn closeness_sq_to_reduction(&self) -> Option<F> {
        None
    }
    fn boxed_split(
        &self,
    ) -> Option<(
        Box<dyn BoxedBezier<F, [F; D]>>,
        Box<dyn BoxedBezier<F, [F; D]>>,
    )> {
        let (b0, b1) = <Self as BezierSplit<_>>::split(self);
        Some((Box::new(b0), Box::new(b1)))
    }
}

impl<F: Num, const D: usize> BezierElevate<F, [F; D]> for [[F; D]; 2] {
    type ElevatedByOne = [[F; D]; 3];
    // Full elevation is not supported
    type Elevated = Self;
    fn elevate_by_one(&self) -> Option<[[F; D]; 3]> {
        Some([
            self[0],
            vector::sum_scaled(self, &[F::frac(1, 2), F::frac(1, 2)]),
            self[1],
        ])
    }
}

impl<F: Num, const D: usize> BezierReduce<F, [F; D]> for [[F; D]; 2] {
    type Reduced = Self;
    type Quadratic = Self;
    type Cubic = Self;
    fn reduce(&self, _method: BezierReduction) -> Option<Self::Reduced> {
        None
    }
    fn can_reduce(&self, _method: BezierReduction) -> bool {
        false
    }
    fn dc_sq_from_reduction(&self, _method: BezierReduction) -> F {
        F::ZERO
    }

    fn dc_sq_from_quadratic(&self) -> F {
        F::ZERO
    }
    fn dc_sq_from_cubic(&self) -> F {
        F::ZERO
    }

    fn reduced_to_quadratic(&self) -> Option<Self::Quadratic> {
        None
    }
    fn reduced_to_cubic(&self) -> Option<Self::Cubic> {
        None
    }
}

impl<F: Num, const D: usize> BezierConstruct<F, D> for [[F; D]; 2] {
    fn of_builder(builder: &BezierBuilder<F, D>) -> Result<Self, BezierError> {
        let mut matrix = [F::ZERO; 4];
        let mut pts = [[F::ZERO; D]; 2];
        builder.fill_fwd_matrix_and_pts(&mut matrix, &mut pts)?;
        if geo_nd::matrix::determinant2(&matrix).is_unreliable_divisor() {
            Err(BezierError::BadBuildConstraints)
        } else {
            let matrix = geo_nd::matrix::inverse2(&matrix);
            let mut bezier = [[F::ZERO; D]; 2];
            bezier[0] = geo_nd::vector::sum_scaled(&pts, &matrix[0..2]);
            bezier[1] = geo_nd::vector::sum_scaled(&pts, &matrix[2..4]);
            Ok(bezier)
        }
    }
}
