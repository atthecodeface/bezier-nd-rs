use crate::utils;
use crate::Num;
use crate::{
    bernstein_fns, BezierDistance, BezierEval, BezierMinMax, BezierReduce, BezierSection,
    BezierSplit, BoxedBezier,
};

use geo_nd::vector;

impl<F: 'static + Num, const D: usize> BezierEval<F, [F; D]> for [[F; D]; 4] {
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
    fn closeness_sq_to_line(&self) -> F {
        self.dc_sq_from_line()
    }
    fn dc_sq_from_line(&self) -> F {
        let one_third: F = (0.33333333).into();
        let two_thirds: F = 0.666_666_7.into();
        let dv_0 = vector::sum_scaled(self, &[-two_thirds, F::ONE, F::ZERO, -one_third]);
        let dc2_0 = vector::length_sq(&dv_0);
        let dv_1 = vector::sum_scaled(self, &[-one_third, F::ZERO, F::ONE, -two_thirds]);
        let dc2_1 = vector::length_sq(&dv_1);
        if dc2_0 < dc2_1 {
            dc2_1
        } else {
            dc2_0
        }
    }
    fn num_control_points(&self) -> usize {
        4
    }
    fn control_point(&self, n: usize) -> &[F; D] {
        &self[n]
    }
}

impl<F: Num, const D: usize> BezierMinMax<F> for [[F; D]; 4] {
    fn t_coord_at_min_max(&self, use_max: bool, pt_index: usize) -> Option<(F, F)> {
        // poly is u^3 P0 + 3u^2t P1 +  3ut^2 P2 + t^3 P3
        //
        // grad/3 is -u^2 P0 + u(1-3t)P1 + t(2-3t).P2 + t^2 P3 = 0
        //
        // i.e. (-1 + 2t - t^2)P0 + (1 + 3t^2 -4t)P1 + (2t - 3t^2)P2 + t^2P3 = 0
        // (-P0 + 3P1 -3P2 + P3)t^2 + 2(P0 -2P1 + P2)t + P1-P0 = 0
        //
        // t = (-(P0 -2P1 + P2) +- sqrt((P0 -2P1 + P2)*(P0 -2P1 + P2)-(P1-P0)(-P0 + 3P1 -3P2 + P3))) / (-P0 + 3P1 -3P2 + P3)
        //
        // If this has real roots then find one with NR; then other is -b/a - first root
        // t = -b/2a +- sqrt()/2a => t0+t1 = -b/a =>
        let p0 = self[0][pt_index];
        let p1 = self[1][pt_index];
        let p2 = self[2][pt_index];
        let p3 = self[3][pt_index];
        let a = p3 - p0 + (p1 - p2) * (3.0_f32).into();
        let p02_sel = utils::min_or_max(use_max, F::ZERO, p0, F::ONE, p3);
        if a.is_unreliable_divisor() {
            todo!();
        }
        let half_b = p0 + p2 - p1 - p1;
        let c = p1 - p0;
        if half_b * half_b < a * c {
            Some(p02_sel)
        } else {
            let b = half_b + half_b;
            let poly = [c, b, a];
            use crate::PolyNewtonRaphson;
            let Some(t0) = poly.find_root_nr(F::ZERO, (1E-6_f32).into()) else {
                return Some(p02_sel);
            };
            let t1 = -(b / a + t0);
            let p02_sel = {
                if t0 > F::ZERO && t0 < F::ONE {
                    let u0 = F::ONE - t0;
                    utils::min_or_max(
                        use_max,
                        p02_sel.0,
                        p02_sel.1,
                        t0,
                        u0 * u0 * u0 * p0
                            + u0 * t0 * (3.0_f32).into() * (u0 * p1 + t0 * p2)
                            + t0 * t0 * t0 * p3,
                    )
                } else {
                    p02_sel
                }
            };
            if t1 > F::ZERO && t1 < F::ONE {
                let u1 = F::ONE - t1;
                Some(utils::min_or_max(
                    use_max,
                    p02_sel.0,
                    p02_sel.1,
                    t1,
                    u1 * u1 * u1 * p0
                        + u1 * t1 * (3.0_f32).into() * (u1 * p1 + t1 * p2)
                        + t1 * t1 * t1 * p3,
                ))
            } else {
                Some(p02_sel)
            }
        }
    }
}

impl<F: Num, const D: usize> BezierDistance<F, [F; D]> for [[F; D]; 4] {
    /// The closest point on a cubic bezier to a point is not analytically determinable
    fn t_dsq_closest_to_pt(&self, _pt: &[F; D]) -> Option<(F, F)> {
        None
    }

    /// An estimate of the minimum distance squared from the Bezier to the point
    /// for 0<=t<=1. This must ALWAYS be less than or equal to the true minimum distance.
    /// If a Bezier does not support this estimate then it *can* return ZERO.
    ///
    /// The points on the Bezier are all within the convex hull of the Bezier, which is
    /// not simple to determine for a cubic.
    ///
    /// The points on the Bezier are also within a distance max(dc_sq) of the line betwen the endpoints
    /// of the Bezier where `dc_sq[i]` is the square of the distance between control point `[i]` and the
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
            // dout <= dmin = sqrt(d_sq) - sqrt(dc_sq)
            // dout^2 <= (sqrt(d_sq) - sqrt(dc_sq)).sqrt(d_sq) - sqrt(dc_sq)
            // dout^2 <= d_sq + dc_sq - 2.sqrt(dc_sq).sqrt(d_sq)
            // dout^2 <= d_sq - dc_sq
            d_sq - dc_sq
        }
    }
}

impl<F: Num, const D: usize> BezierSplit for [[F; D]; 4] {
    fn split(&self) -> (Self, Self) {
        let pm = vector::sum_scaled(
            self,
            &[
                (0.125_f32).into(),
                (0.375_f32).into(),
                (0.375_f32).into(),
                (0.125_f32).into(),
            ],
        );
        let c00 = vector::sum_scaled(
            self,
            &[
                (0.5_f32).into(),
                (0.5_f32).into(),
                (0.0_f32).into(),
                (0.0_f32).into(),
            ],
        );
        let c01 = vector::sum_scaled(
            self,
            &[
                (0.25_f32).into(),
                (0.5_f32).into(),
                (0.25_f32).into(),
                (0.0_f32).into(),
            ],
        );
        let c10 = vector::sum_scaled(
            self,
            &[
                (0.0_f32).into(),
                (0.25_f32).into(),
                (0.5_f32).into(),
                (0.25_f32).into(),
            ],
        );
        let c11 = vector::sum_scaled(
            self,
            &[
                (0.0_f32).into(),
                (0.0_f32).into(),
                (0.5_f32).into(),
                (0.5_f32).into(),
            ],
        );
        ([self[0], c00, c01, pm], [pm, c10, c11, self[3]])
    }
}

impl<F: Num, const D: usize> BezierSection<F> for [[F; D]; 4] {
    /// Split the Bezier into two (one covering parameter values 0..t, the other t..1)
    fn split_at(&self, t: F) -> (Self, Self) {
        let mut latter = *self;
        let mut first = [[F::ZERO; D]; 4];
        bernstein_fns::split::into_two_at_de_cast(&mut latter, t, &mut first);
        (first, latter)
    }

    /// Generate the Bezier with parameter t' where 0<=t'<=1 provides the same points
    /// as the original Bezier with t0<=t<=1
    fn section(&self, t0: F, t1: F) -> Self {
        let mut to_split = *self;
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

impl<F: Num, const D: usize> BoxedBezier<F, [F; D]> for [[F; D]; 4] {
    fn closeness_sq_to_reduction(&self) -> Option<F> {
        None
    }
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

impl<F: 'static + Num, const D: usize> BezierReduce<F, [F; D]> for [[F; D]; 4] {
    type Reduced = [[F; D]; 3];
    type Quadratic = [[F; D]; 3];
    type Cubic = Self;
    fn reduce(&self) -> Self::Reduced {
        todo!();
    }
    fn can_reduce(&self) -> bool {
        true
    }
    fn closeness_sq_to_reduction(&self) -> Option<F> {
        None
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

    fn reduced_to_quadratic(&self) -> Option<Self::Quadratic> {
        todo!()
    }
    fn reduced_to_cubic(&self) -> Option<Self::Cubic> {
        None
    }
}
