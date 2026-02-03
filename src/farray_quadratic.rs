use crate::utils;
use crate::{
    bernstein_fns, BezierBuilder, BezierConstruct, BezierDistance, BezierEval, BezierMinMax,
    BezierReduce, BezierSection, BezierSplit, BoxedBezier,
};
use crate::{Float, Num};
use geo_nd::vector;

impl<F: Num, const D: usize> BezierEval<F, [F; D]> for [[F; D]; 3] {
    fn point_at(&self, t: F) -> [F; D] {
        let two: F = (2.0_f32).into();
        let u = F::ONE - t;
        vector::sum_scaled(self, &[u * u, two * u * t, t * t])
    }
    fn derivative_at(&self, t: F) -> (F, [F; D]) {
        let two: F = (2.0_f32).into();
        let u = F::ONE - t;
        (two, vector::sum_scaled(self, &[-u, u - t, t]))
    }
    fn endpoints(&self) -> (&[F; D], &[F; D]) {
        (&self[0], &self[2])
    }
    fn closeness_sq_to_line(&self) -> F {
        eprintln!("quad {self:?}");
        utils::distance_sq_to_line_segment(&self[1], &self[0], &self[2])
    }
    fn dc_sq_from_line(&self) -> F {
        eprintln!("dc {self:?}");
        vector::length_sq(&vector::sum_scaled(
            self,
            &[(0.5_f32).into(), -F::ONE, (0.5_f32).into()],
        ))
    }

    fn num_control_points(&self) -> usize {
        3
    }
    fn control_point(&self, n: usize) -> &[F; D] {
        &self[n]
    }
}

impl<F: Num, const D: usize> BezierMinMax<F> for [[F; D]; 3] {
    fn t_coord_at_min_max(&self, use_max: bool, pt_index: usize) -> Option<(F, F)> {
        // poly is u^2 P0 + 2ut P1 + t^2 P2
        //
        // grad/2 is (t-1) P0 + (1-2t) P1 + t P2 = 0 if t(P0-2P1+P2) = P0-P1
        //
        // i.e. t = (P0-P1) / (P0-2P1+P2)
        //
        // Second differential is P0 -2P1 + P2
        let p0 = self[0][pt_index];
        let p1 = self[1][pt_index];
        let p2 = self[2][pt_index];
        let p02_sel = utils::min_or_max(use_max, F::ZERO, p0, F::ONE, p2);
        let d_dt_denom = p0 + p2 - p1 * (2.0_f32).into();
        let d_dt_numer = p0 - p1;
        if d_dt_denom.is_unreliable_divisor() {
            Some(p02_sel)
        } else {
            let t = d_dt_numer / d_dt_denom;
            if t > F::ZERO && t < F::ONE {
                let u = F::ONE - t;
                Some(utils::min_or_max(
                    use_max,
                    p02_sel.0,
                    p02_sel.1,
                    t,
                    u * u * p0 + u * t * p1 * (2.0_f32).into() + t * t * p2,
                ))
            } else {
                Some(p02_sel)
            }
        }
    }
}

impl<F: Float, const D: usize> BezierDistance<F, [F; D]> for [[F; D]; 3] {
    /// A value of t and distance squared to it, 0<=t<=1, for which the distance between the point
    /// P and the quadratic Bezier between the two points is a minimum.
    ///
    /// The distance squared between a quadratic Bezier and a point P is a quartic equation in t (with
    /// a positive coefficient of t^4)
    /// and the derivative is a cubic which has at least one real root (as there is a least one minimum
    /// in the distance squared); indeed, if there are three real roots then they can be ordered t0<t1<t2,
    /// and the minima are at t0 and t2.
    ///
    /// This calculates the roots, the distances to the endpoints and the minima roots, and returns the smallest
    /// distance squared and t value.
    // Note p(t) = u^2.p0 + 2ut.p1 + t^2.p2
    // or p0 + 2(p1 - p0)t + (p0 - 2p1 + p2).t^2
    // and
    // p'(t) = 2(p1-p0) + 2(p0-2p1+p2)t
    //
    // Define p012 = p0-2*p1+p2, p01=p1-p0
    //
    // dsq = distance_sq to point = |p(t)-p|^2 = (p(t)-p) . (p(t)-p)
    // Define p0p = p0-p
    //
    // d/dt (dsq) = 2 (p(t)-p) . d/dt(p(t)-p) = 2 (p(t)-p) . p'(t)
    // = 0 at a minimum distance for the curve to p
    //
    // (p(t)-p) . p'(t) = (p0 + 2(p1 - p0)t + (p0 - 2p1 + p2).t^2 - p) . (2(p1-p0) + 2(p0-2p1+p2)t)
    // = (p0p + 2p01.t + p012.t^2) . 2 . (p01 + p012.t)
    // = 2(p012.p012.t^3 + (2p01.p012 + p012.p01)t^2 + (2.p01.p01 + p0p.p012)t + p0p.p01)
    // = 2(p012.p012.t^3 + 3.p01.p012.t^2 + (2.p01.p01 + p0p.p012)t + p0p.p01)
    // = 0 at minimum or maximum distance squared for varying t
    //
    // This cubic can be solved for t, leading to one or three roots
    //
    // d_sq is a quartic that is always positive; it will have a positive coefficient of t^4. It
    // s derivative will have either one or three real roots; if there are three then they can be
    // sorted (into t0, t1, t2; t0 <= t1 <= t2 and the minima of d_sq will be at t0 and t2 (but
    // the closest points on 0<=t<=1 may be at t=0 or t=1)
    //
    fn t_dsq_closest_to_pt(&self, pt: &[F; D]) -> Option<(F, F)> {
        // First find the closest endpoint in case we have to bail (rather than panic)
        let d0_sq = vector::distance_sq(pt, &self[0]);
        let d1_sq = vector::distance_sq(pt, &self[2]);
        let mut t_min = F::ZERO;
        let mut dsq_min = d0_sq;
        if d1_sq < dsq_min {
            t_min = F::ONE;
            dsq_min = d1_sq;
        }

        let p012 = vector::add(
            vector::add(self[0], &self[1], (-2.0_f32).into()),
            &self[2],
            F::ONE,
        );
        let p01 = vector::add(self[1], &self[0], (-1.0_f32).into());
        let p0p = vector::add(self[0], pt, (-2.0_f32).into());
        let mut d_dt_dsq_poly = [
            vector::dot(&p01, &p0p),
            vector::dot(&p01, &p01) * (2.0_f32).into() + vector::dot(&p012, &p0p),
            vector::dot(&p01, &p012) * (3.0_f32).into(),
            vector::dot(&p012, &p012),
        ];
        use crate::polynomial::{PolyFindRoots, PolyNewtonRaphson, Polynomial};
        let Some(t0) = d_dt_dsq_poly.find_root_nr(0.5_f32.into(), 1E7_f32.into()) else {
            // panic!("Failed to find root for polynomial {d_dt_dsq_poly:?}");
            return Some((t_min, dsq_min));
        };
        let mut d_dt_dsq_poly_without_root = [F::ZERO; 3];
        if !d_dt_dsq_poly_without_root.set_divide(
            &mut d_dt_dsq_poly,
            &[-t0, F::ONE],
            (1.0E-4_f32).into(),
        ) {
            // panic!("Failed to divide by x-root!");
            return Some((t_min, dsq_min));
        }
        let (opt_t1, opt_t2) = d_dt_dsq_poly_without_root.find_roots_quad();
        let t1 = opt_t1.unwrap_or(t0);
        let t2 = opt_t2.unwrap_or(t0);
        // Make t0 < t1; t1 and t2 will be unordered, and t0 and t2 will be unordered
        let (t0, t1) = {
            if t0 < t1 {
                (t0, t1)
            } else {
                (t1, t0)
            }
        };
        // Make t0 < t2 ordered; t0 < t1 and t0<t2; ordering of t1 and t2 uncertain
        let (t0, t2) = {
            if t0 < t2 {
                (t0, t2)
            } else {
                (t2, t0)
            }
        };
        // Make t0 < t2 ordered; t0 < t1 and t0<t2; ordering of t1 and t2 uncertain
        let (t1, t2) = {
            if t1 < t2 {
                (t1, t2)
            } else {
                (t2, t1)
            }
        };
        // t1 is not a minimum, so we don't use that any further
        let _t1 = t1;

        // At this point t0 and t2 are minima, and t0<t2
        //
        // There are six cases:
        //         0.             1.           Minimum
        // --------------------------------------------
        //  t0 t2  |              |          | at t=0
        //  t0     |  t2          |          | at t=0 or t2
        //  t0     |              |  t2      | at t=0 or t=1
        //         |  t0 t2       |          | at t0 or t2
        //         |  t0          |  t2      | at t0 or t=1
        //         |              |  t0 t2.  } at t=1
        let t0_lt_zero = t0 <= F::ZERO;
        let t2_lt_zero = t2 <= F::ZERO;
        let t0_gt_one = t0 >= F::ONE;
        let t2_gt_one = t2 >= F::ONE;
        if t2_lt_zero || t0_gt_one || (t0_lt_zero && t2_gt_one) {
            return Some((t_min, dsq_min));
        }
        let dsq_at_t0 = vector::distance_sq(pt, &self.point_at(t0));
        let dsq_at_t1 = vector::distance_sq(pt, &self.point_at(t1));
        if dsq_at_t0 < dsq_min {
            dsq_min = dsq_at_t0;
            t_min = t0;
        }
        if dsq_at_t1 < dsq_min {
            dsq_min = dsq_at_t1;
            t_min = t1;
        }
        Some((t_min, dsq_min))
    }

    /// An estimate of the minimum distance squared from the Bezier to the point
    /// for 0<=t<=1.
    ///
    /// This returns the minimum distance betwen the point and the convex hull of the Bezier,
    /// which is the triangle of the three control points
    ///
    /// If the point projected onto the triangle plane is within the triangle then this is the
    /// distance squared between the projected point and the point
    ///
    /// If the point projected onto the triangle plane is outside the triangle then it cannot
    /// be closer to the Bezier than it is from the closest edge of the triangle, so the minimum
    /// distance squared between the Point and the three triangle edges will suffice
    fn est_min_distance_sq_to(&self, pt: &[F; D]) -> F {
        if let Some((k1, k2)) =
            crate::utils::barycentric_coordinates(pt, &self[0], &self[1], &self[2])
        {
            let k0 = F::ONE - k1 - k2;
            if k0 >= F::ZERO
                && k1 >= F::ZERO
                && k2 >= F::ZERO
                && k0 <= F::ONE
                && k1 <= F::ONE
                && k2 <= F::ONE
            {
                // For 2-dimensional Beziers the point projected *is* the point
                if D <= 2 {
                    return F::ZERO;
                } else {
                    let pt_rel = vector::sub(*pt, &self[0], k0);
                    let pt_rel = vector::sub(pt_rel, &self[1], k1);
                    let pt_rel = vector::sub(pt_rel, &self[2], k2);
                    return vector::length_sq(&pt_rel);
                }
            }
        }
        let d0 = crate::utils::distance_sq_to_line_segment(pt, &self[0], &self[1]);
        let d1 = crate::utils::distance_sq_to_line_segment(pt, &self[0], &self[2]);
        let d2 = crate::utils::distance_sq_to_line_segment(pt, &self[1], &self[2]);
        d0.min(d1.min(d2))
    }
}

impl<F: Num, const D: usize> BezierSplit for [[F; D]; 3] {
    fn split(&self) -> (Self, Self) {
        let c0 = vector::sum_scaled(self, &[0.5_f32.into(), 0.5_f32.into(), F::ZERO]);
        let c1 = vector::sum_scaled(self, &[F::ZERO, 0.5_f32.into(), 0.5_f32.into()]);
        let pm = vector::sum_scaled(self, &[0.25_f32.into(), 0.5_f32.into(), 0.25_f32.into()]);
        ([self[0], c0, pm], [pm, c1, self[2]])
    }
}

impl<F: Num, const D: usize> BezierSection<F> for [[F; D]; 3] {
    /// Split the Bezier into two (one covering parameter values 0..t, the other t..1)
    fn split_at(&self, t: F) -> (Self, Self) {
        let mut latter = *self;
        let mut first = [[F::ZERO; D]; 3];
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

// This requires Float as BezierDistance for [[F;D];3] requires Float for closest point to curve
impl<F: Num, const D: usize> BoxedBezier<F, [F; D]> for [[F; D]; 3] {
    fn boxed_reduce(&self) -> Option<Box<dyn BoxedBezier<F, [F; D]>>> {
        Some(Box::new([self[0], self[1]]))
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
        let (b0, b1) = <Self as BezierSplit>::split(self);
        Some((Box::new(b0), Box::new(b1)))
    }
}

impl<F: Num, const D: usize> BezierReduce<F, [F; D]> for [[F; D]; 3] {
    type Reduced = [[F; D]; 2];
    type Quadratic = Self;
    type Cubic = Self;
    fn reduce(&self) -> Self::Reduced {
        [self[0], self[1]]
    }
    fn can_reduce(&self) -> bool {
        true
    }
    fn closeness_sq_to_reduction(&self) -> Option<F> {
        todo!()
    }

    fn closeness_sq_to_quadratic(&self) -> F {
        F::ZERO
    }
    fn closeness_sq_to_cubic(&self) -> F {
        F::ZERO
    }
    fn reduced_to_quadratic(&self) -> Option<Self::Quadratic> {
        None
    }
    fn reduced_to_cubic(&self) -> Option<Self::Cubic> {
        None
    }
}

impl<F: Num, const D: usize> BezierConstruct<F, D> for [[F; D]; 3] {
    fn of_builder(builder: &BezierBuilder<F, D>) -> Result<Self, ()> {
        let mut matrix = [F::ZERO; 9];
        let mut pts = [[F::ZERO; D]; 3];
        builder.fill_fwd_matrix_and_pts(&mut matrix, &mut pts)?;
        if geo_nd::matrix::determinant3(&matrix).is_unreliable_divisor() {
            Err(())
        } else {
            let matrix = geo_nd::matrix::inverse3(&matrix);
            let mut bezier = [[F::ZERO; D]; 3];
            bezier[0] = geo_nd::vector::sum_scaled(&pts, &matrix[0..3]);
            bezier[1] = geo_nd::vector::sum_scaled(&pts, &matrix[3..6]);
            bezier[2] = geo_nd::vector::sum_scaled(&pts, &matrix[6..9]);
            Ok(bezier)
        }
    }
}
