use crate::utils;
use crate::{
    bernstein_fns, BezierBuilder, BezierConstruct, BezierDistance, BezierElevate, BezierEval,
    BezierMinMax, BezierOps, BezierReduce, BezierSection, BezierSplit, BoxedBezier,
};

use crate::Num;
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
        utils::distance_sq_to_line_segment(&self[1], &self[0], &self[2])
    }
    fn dc_sq_from_line(&self) -> F {
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
    fn for_each_control_point(&self, map: &mut dyn FnMut(usize, &[F; D])) {
        self.iter().enumerate().for_each(|(i, pt)| map(i, pt))
    }
}

impl<F: Num, const D: usize> BezierOps<F, [F; D]> for [[F; D]; 3] {
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

impl<F: Num, const D: usize> BezierDistance<F, [F; D]> for [[F; D]; 3] {
    fn t_dsq_closest_to_pt(&self, pt: &[F; D]) -> Option<(F, F)> {
        bernstein_fns::distance::bezier_quad_t_dsq_closest_to_pt(self, pt)
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
        crate::utils::min(d0, crate::utils::min(d1, d2))
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
        bernstein_fns::de_casteljau::split_at(&mut latter, t, &mut first);
        (first, latter)
    }

    /// Generate the Bezier with parameter t' where 0<=t'<=1 provides the same points
    /// as the original Bezier with t0<=t<=1
    fn section(&self, t0: F, t1: F) -> Self {
        let mut to_split = *self;
        if t0 > F::ZERO {
            bernstein_fns::de_casteljau::bezier_from(&mut to_split, t0);
        }
        if t1 < F::ONE {
            let t10 = (t1 - t0) / (F::ONE - t0);
            bernstein_fns::de_casteljau::bezier_to(&mut to_split, t10);
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

impl<F: Num, const D: usize> BezierElevate<F, [F; D]> for [[F; D]; 3] {
    type ElevatedByOne = [[F; D]; 4];
    // Full elevation is not supported
    type Elevated = Self;
    fn elevate_by_one(&self) -> Option<[[F; D]; 4]> {
        Some([
            self[0],
            vector::sum_scaled(&self[0..2], &[0.3333333_f32.into(), 0.66666667_f32.into()]),
            vector::sum_scaled(&self[1..3], &[0.66666667_f32.into(), 0.3333333_f32.into()]),
            self[2],
        ])
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
