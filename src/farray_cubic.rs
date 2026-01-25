use geo_nd::{vector, Float, Num};

use crate::{BezierDistance, BezierEval, BezierReduce, BezierSplit, BoxedBezier};

impl<F: 'static + Num + From<f32>, const D: usize> BezierEval<F, [F; D]> for [[F; D]; 4] {
    fn point_at(&self, t: F) -> [F; D] {
        let three: F = (3.0_f32).into();
        let u = F::ONE - t;
        vector::sum_scaled(
            self,
            &[u * u * u, three * u * u * t, three * u * t * t, t * t * t],
        )
    }
    fn derivative_at(&self, t: F) -> (F, [F; D]) {
        let three: F = (3.0_f32).into();
        let u = F::ONE - t;
        (
            three,
            vector::sum_scaled(self, &[-u * u, u * (u - t - t), t * (u + u - t), t * t]),
        )
    }
    fn endpoints(&self) -> (&[F; D], &[F; D]) {
        (&self[0], &self[3])
    }
    fn is_straight(&self, straightness_sq: F) -> bool {
        let dc_sq = crate::utils::straightness_sq_of_cubic(self);
        dc_sq < straightness_sq
    }
    fn closeness_sq_to_quadratic(&self) -> F {
        let m_half = (-0.5_f32).into();
        let dv_0 = vector::sum_scaled(self, &[m_half, F::ONE, F::ZERO, m_half]);
        let dc2_0 = vector::length_sq(&dv_0);
        let dv_1 = vector::sum_scaled(self, &[m_half, F::ZERO, F::ONE, m_half]);
        let dc2_1 = vector::length_sq(&dv_1);
        if dc2_0 < dc2_1 {
            dc2_1
        } else {
            dc2_0
        }
    }
    fn closeness_sq_to_cubic(&self) -> F {
        F::ZERO
    }
    fn num_control_points(&self) -> usize {
        4
    }
    fn control_point(&self, n: usize) -> &[F; D] {
        &self[n]
    }
}

impl<F: 'static + Float, const D: usize> BezierDistance<F, [F; D]> for [[F; D]; 4] {
    /// The closest point on a cubic bezier to a point is not analytically determinable
    fn t_dsq_closest_to_pt(&self, pt: &[F; D]) -> Option<(F, F)> {
        None
    }

    /// An estimate of the minimum distance squared from the Bezier to the point
    /// for 0<=t<=1.
    ///
    /// The points on the Bezier are all within the convex hull of the Bezier, which is
    /// not simple to determine for a cubic.
    ///
    /// The points on the Bezier are also within a distance max(dc_sq) of the line betwen the endpoints
    /// of the Bezier where dc_sq[i] is the square of the distance between control point [i] and the
    /// position of that control point on a linear Bezier elevated to the degree of this Bezier.
    ///
    /// Hence if the point is at a distance squared of d_sq from the linear Bezier, then it is potentially
    /// at a distance of d_sq-max(dc_sq) (or 0 if that is negative) away from the Bezier curve.
    fn est_min_distance_sq_to(&self, pt: &[F; D]) -> F {
        let d_sq = crate::utils::distance_sq_to_line_segment(pt, &self[0], &self[3]);
        let dc_sq = crate::utils::straightness_sq_of_cubic(self);
        if d_sq < dc_sq {
            F::ZERO
        } else {
            let dc = dc_sq.sqrt();
            let d = d_sq.sqrt();
            (dc - d) * (dc - d)
        }
    }
}

impl<F: 'static + Num + From<f32>, const D: usize> BezierSplit for [[F; D]; 4] {
    fn split(&self) -> (Self, Self) {
        todo!();
    }
}

impl<F: 'static + Num + From<f32>, const D: usize> BoxedBezier<F, [F; D]> for [[F; D]; 4] {
    fn boxed_reduce(&self) -> Option<Box<dyn BoxedBezier<F, [F; D]>>> {
        Some(Box::new([self[0], self[1]]))
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

impl<F: 'static + Num + From<f32>, const D: usize> BezierReduce<F, [F; D]> for [[F; D]; 4] {
    type Reduced = [[F; D]; 3];
    type Quadratic = [[F; D]; 3];
    type Cubic = Self;
    fn reduce(&self) -> Self::Reduced {
        todo!();
    }
    fn can_reduce() -> bool {
        true
    }
    fn closeness_sq_to_reduction(&self) -> F {
        todo!()
    }

    fn reduced_to_quadratic(&self) -> Option<Self::Quadratic> {
        todo!()
    }
    fn reduced_to_cubic(&self) -> Option<Self::Cubic> {
        None
    }
}
