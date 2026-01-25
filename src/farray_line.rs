use geo_nd::{vector, Num};

use crate::{BezierDistance, BezierEval, BezierSplit, BoxedBezier};

impl<F: 'static + Num + From<f32>, const D: usize> BezierEval<F, [F; D]> for [[F; D]; 2] {
    fn point_at(&self, t: F) -> [F; D] {
        vector::add(vector::scale(self[0], F::ONE - t), &self[1], t)
    }
    fn derivative_at(&self, _t: F) -> (F, [F; D]) {
        (F::ONE, vector::add(self[1], &self[0], -F::ONE))
    }
    fn endpoints(&self) -> (&[F; D], &[F; D]) {
        (&self[0], &self[1])
    }
    fn is_straight(&self, _straightness_sq: F) -> bool {
        true
    }
    fn closeness_sq_to_quadratic(&self) -> F {
        F::ZERO
    }
    fn closeness_sq_to_cubic(&self) -> F {
        F::ZERO
    }
    fn num_control_points(&self) -> usize {
        2
    }
    fn control_point(&self, n: usize) -> &[F; D] {
        &self[n]
    }
}

impl<F: 'static + Num + From<f32>, const D: usize> BezierSplit for [[F; D]; 2] {
    fn split(&self) -> (Self, Self) {
        let m = vector::scale(
            vector::add(self[0], &self[1], 1.0_f32.into()),
            0.5_f32.into(),
        );
        ([self[0], m], [m, self[1]])
    }
}

impl<F: 'static + Num + From<f32>, const D: usize> BoxedBezier<F, [F; D]> for [[F; D]; 2] {
    fn boxed_reduce(&self) -> Option<Box<dyn BoxedBezier<F, [F; D]>>> {
        None
    }
    fn boxed_split(
        &self,
    ) -> Option<(
        Box<dyn BoxedBezier<F, [F; D]>>,
        Box<dyn BoxedBezier<F, [F; D]>>,
    )> {
        let (b0, b1) = <Self as BezierSplit>::split(self);
        Some((Box::new(b0), Box::new(b1)))
    }
}

impl<F: 'static + geo_nd::Float, const D: usize> BezierDistance<F, [F; D]> for [[F; D]; 2] {
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
        let (t_times_len_sq, line_len_sq, valid) =
            crate::utils::relative_to_line(pt, &self[0], &self[1]);
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
        crate::utils::distance_sq_to_line_segment(p, &self[0], &self[1])
    }
}
