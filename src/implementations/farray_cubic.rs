use crate::utils;
use crate::Num;
use crate::{
    bernstein_fns, BezierBuilder, BezierConstruct, BezierElevate, BezierEval, BezierFlatIterator,
    BezierLineIter, BezierLineTIter, BezierOps, BezierReduce, BezierReduction, BezierSplit,
    BoxedBezier,
};

use geo_nd::vector;

impl<F: Num, const D: usize> BezierEval<F, [F; D]> for [[F; D]; 4] {
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
    fn endpoints(&self) -> ([F; D], [F; D]) {
        (self[0], self[3])
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
    fn control_points(&self) -> &[[F; D]] {
        self
    }
    fn t_coords_at_min_max(
        &self,
        pt_index: usize,
        give_min: bool,
        give_max: bool,
    ) -> (Option<(F, F)>, Option<(F, F)>) {
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
        let (mut opt_min, mut opt_max) =
            utils::opt_min_and_max_tc(give_min, give_max, (F::ZERO, p0), (F::ONE, p3), None);
        let a = p3 - p0 + (p1 - p2) * (3.0_f32).into();
        let b = (p0 + p2 - p1 - p1) * 2.0_f32.into();
        let c = p1 - p0;
        let (opt_t0, opt_t1) = utils::find_real_roots_quad_num(&[c, b, a]);
        if let Some(t0) = opt_t0 {
            if t0 > F::ZERO && t0 < F::ONE {
                let c0 = self.point_at(t0)[pt_index];
                opt_min = opt_min.map(|t_c| utils::min_tc(t_c, (t0, c0)));
                opt_max = opt_max.map(|t_c| utils::max_tc(t_c, (t0, c0)));
            }
        }
        if let Some(t1) = opt_t1 {
            if t1 > F::ZERO && t1 < F::ONE {
                let c1 = self.point_at(t1)[pt_index];
                opt_min = opt_min.map(|t_c| utils::min_tc(t_c, (t1, c1)));
                opt_max = opt_max.map(|t_c| utils::max_tc(t_c, (t1, c1)));
            }
        }
        (opt_min, opt_max)
    }

    /// The closest point on a cubic bezier to a point is not analytically determinable
    fn t_dsq_closest_to_pt(&self, _pt: &[F; D]) -> Option<(F, F)> {
        None
    }

    fn est_min_distance_sq_to(&self, pt: &[F; D]) -> F {
        let d_sq = utils::distance_sq_to_line_segment(pt, &self[0], &self[3]);
        let dc_sq = utils::straightness_sq_of_cubic(self);
        if d_sq < dc_sq {
            F::ZERO
        } else {
            utils::est_d_m_c_from_dsq_m_dcsq(d_sq, dc_sq)
        }
    }
}

impl<F: Num, const D: usize> BezierOps<F, [F; D]> for [[F; D]; 4] {
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
    fn map_all_pts(&mut self, map: &dyn Fn(&mut [[F; D]]) -> bool) -> bool {
        map(self)
    }
}

impl<F: Num, const D: usize> BezierSplit<F> for [[F; D]; 4] {
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

    /// Split the Bezier into two (one covering parameter values 0..t, the other t..1)
    fn split_at(&self, t: F) -> (Self, Self) {
        let mut latter = *self;
        let mut first = [[F::ZERO; D]; 4];
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

impl<F, const D: usize> BezierFlatIterator<F, [F; D]> for [[F; D]; 4]
where
    F: Num,
{
    fn as_lines(&self, closeness_sq: F) -> impl Iterator<Item = ([F; D], [F; D])> {
        BezierLineIter::<_, _, _, false>::new(self, closeness_sq)
    }
    fn as_t_lines(&self, closeness_sq: F) -> impl Iterator<Item = (F, [F; D], F, [F; D])> {
        BezierLineTIter::<_, _, _, false>::new(self, closeness_sq)
    }
    fn as_t_lines_dc(&self, closeness_sq: F) -> impl Iterator<Item = (F, [F; D], F, [F; D])> {
        BezierLineTIter::<_, _, _, true>::new(self, closeness_sq)
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
        let (b0, b1) = <Self as BezierSplit<_>>::split(self);
        Some((Box::new(b0), Box::new(b1)))
    }
}

impl<F: Num, const D: usize> BezierElevate<F, [F; D]> for [[F; D]; 4] {
    type ElevatedByOne = Self;
    type Elevated = Self;
    fn elevate_by_one(&self) -> Option<[[F; D]; 4]> {
        None
    }
}

impl<F: Num, const D: usize> BezierReduce<F, [F; D]> for [[F; D]; 4] {
    type Reduced = [[F; D]; 3];
    type Quadratic = [[F; D]; 3];
    type Cubic = Self;
    fn reduce(&self, method: BezierReduction) -> Self::Reduced {
        todo!();
    }
    fn can_reduce(&self, method: BezierReduction) -> bool {
        true
    }
    fn closeness_sq_to_reduction(&self, method: BezierReduction) -> Option<F> {
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

impl<F: Num, const D: usize> BezierConstruct<F, D> for [[F; D]; 4] {
    fn of_builder(builder: &BezierBuilder<F, D>) -> Result<Self, ()> {
        let mut matrix = [F::ZERO; 16];
        let mut pts = [[F::ZERO; D]; 4];
        builder.fill_fwd_matrix_and_pts(&mut matrix, &mut pts)?;
        if geo_nd::matrix::determinant4(&matrix).is_unreliable_divisor() {
            Err(())
        } else {
            let matrix = geo_nd::matrix::inverse4(&matrix);
            let mut bezier = [[F::ZERO; D]; 4];
            bezier[0] = geo_nd::vector::sum_scaled(&pts, &matrix[0..4]);
            bezier[1] = geo_nd::vector::sum_scaled(&pts, &matrix[4..8]);
            bezier[2] = geo_nd::vector::sum_scaled(&pts, &matrix[8..12]);
            bezier[3] = geo_nd::vector::sum_scaled(&pts, &matrix[12..16]);
            Ok(bezier)
        }
    }
}
