use crate::Num;
use crate::{
    bernstein_fns, constants, metrics, utils, BezierBuilder, BezierConstruct, BezierElevate,
    BezierEval, BezierFlatIterator, BezierLineIter, BezierLineTIter, BezierMetric, BezierOps,
    BezierReduce, BezierReduction, BezierSplit,
};

use geo_nd::{matrix, vector};

impl<F: Num, const D: usize> BezierEval<F, [F; D]> for Vec<[F; D]> {
    fn point_at(&self, t: F) -> [F; D] {
        bernstein_fns::values::point_at(self, t)
    }
    fn derivative_at(&self, t: F) -> (F, [F; D]) {
        bernstein_fns::values::derivative_at(self, t)
    }
    fn closeness_sq_to_line(&self) -> F {
        metrics::dc_sq_from_line(self)
    }
    fn dc_sq_from_line(&self) -> F {
        metrics::dc_sq_from_line(self)
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

    fn t_coords_at_min_max(
        &self,
        pt_index: usize,
        give_min: bool,
        give_max: bool,
    ) -> (Option<(F, F)>, Option<(F, F)>) {
        match self.len() {
            1 => (
                give_min.then_some((F::ZERO, self[0][pt_index])),
                give_max.then_some((F::ZERO, self[0][pt_index])),
            ),
            2 => <&[[F; D]] as TryInto<&[[F; D]; 2]>>::try_into(&self[0..2])
                .unwrap()
                .t_coords_at_min_max(pt_index, give_min, give_max),
            3 => <&[[F; D]] as TryInto<&[[F; D]; 3]>>::try_into(&self[0..3])
                .unwrap()
                .t_coords_at_min_max(pt_index, give_min, give_max),
            4 => <&[[F; D]] as TryInto<&[[F; D]; 4]>>::try_into(&self[0..4])
                .unwrap()
                .t_coords_at_min_max(pt_index, give_min, give_max),
            _ => (None, None),
        }
    }

    fn t_dsq_closest_to_pt(&self, pt: &[F; D]) -> Option<(F, F)> {
        metrics::t_dsq_closest_to_pt(&self, pt)
    }

    fn est_min_distance_sq_to(&self, pt: &[F; D]) -> F {
        if self.len() == 0 {
            F::ZERO
        } else {
            metrics::est_min_distance_sq_to(self, pt)
        }
    }
}

impl<F: Num, const D: usize> BezierOps<F, [F; D]> for Vec<[F; D]> {
    fn add(&mut self, other: &Self) -> bool {
        if self.len() != other.len() {
            false
        } else {
            for (s, o) in self.iter_mut().zip(other.iter()) {
                *s = vector::add(*s, o, F::ONE)
            }
            true
        }
    }

    fn sub(&mut self, other: &Self) -> bool {
        if self.len() != other.len() {
            false
        } else {
            for (s, o) in self.iter_mut().zip(other.iter()) {
                *s = vector::sub(*s, o, F::ONE)
            }
            true
        }
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

impl<F: Num, const D: usize> BezierSplit<F> for Vec<[F; D]> {
    fn split(&self) -> (Self, Self) {
        self.split_at(0.5_f32.into())
    }
    /// Split the Bezier into two (one covering parameter values 0..t, the other t..1)
    fn split_at(&self, t: F) -> (Self, Self) {
        let mut latter = self.clone();
        let mut first = self.clone();
        bernstein_fns::de_casteljau::split_at(&mut latter, t, &mut first);
        (first, latter)
    }

    /// Generate the Bezier with parameter t' where 0<=t'<=1 provides the same points
    /// as the original Bezier with t0<=t<=1
    fn section(&self, t0: F, t1: F) -> Self {
        let mut to_split = self.clone();
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

impl<F, const D: usize> BezierFlatIterator<F, [F; D]> for Vec<[F; D]>
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

impl<F: Num, const D: usize> BezierElevate<F, [F; D]> for Vec<[F; D]> {
    type ElevatedByOne = Self;
    type Elevated = Self;
    fn elevate_by_one(&self) -> Option<Self> {
        let n = self.len();
        let mut s = vec![[F::ZERO; D]; n + 1];
        s[0] = *self.first().unwrap();
        s[n] = *self.last().unwrap();
        let n_f: F = (n as f32).into();
        let mut s0 = F::ONE;
        for (p, (p0, p1)) in s
            .iter_mut()
            .skip(1)
            .zip(self.iter().zip(self.iter().skip(1)))
        {
            *p = vector::sum_scaled(&[*p0, *p1], &[s0 / n_f, F::ONE - s0 / n_f]);
            s0 += F::ONE;
        }
        Some(s)
    }
    fn elevate_by(&self, degree: usize) -> Option<Self> {
        if degree == 0 {
            Some(self.clone())
        } else {
            let mut bezier = self.elevate_by_one();
            for _ in 1..degree {
                bezier = bezier.unwrap().elevate_by_one();
            }
            bezier
        }
    }
}

impl<F, const D: usize> BezierReduce<F, [F; D]> for Vec<[F; D]>
where
    F: Num,
{
    type Reduced = Self;
    type Quadratic = Self;
    type Cubic = Self;

    fn can_reduce(&self, method: BezierReduction) -> bool {
        let table = constants::reduce_table_of_match(method);
        self.len() > 2 && self.len() <= table.len() + 3
    }

    fn reduce(&self, method: BezierReduction) -> Option<Self> {
        match self.len() {
            3 => <&[[F; D]] as TryInto<&[[F; D]; 3]>>::try_into(&self[0..3])
                .unwrap()
                .reduce(method)
                .map(|a| a.into()),
            4 => <&[[F; D]] as TryInto<&[[F; D]; 4]>>::try_into(&self[0..4])
                .unwrap()
                .reduce(method)
                .map(|a| a.into()),
            _ => {
                let table = constants::reduce_table_of_match(method);
                if self.len() > table.len() + 3 {
                    None
                } else {
                    let mut result = vec![[F::ZERO; D]; self.len() - 1];
                    crate::lazy_constants::use_constants_table(
                        |table| bernstein_fns::transform::transform_pts(table, self, &mut result),
                        table,
                        self.len() - 3,
                    );
                    Some(result)
                }
            }
        }
    }

    fn dc_sq_from_reduction(&self, method: BezierReduction) -> F {
        match self.len() {
            3 => <&[[F; D]] as TryInto<&[[F; D]; 3]>>::try_into(&self[0..3])
                .unwrap()
                .dc_sq_from_reduction(method),
            4 => self.dc_sq_from_quadratic(),
            _ => {
                let table = constants::er_minus_i_table_of_match(method);
                if self.len() > table.len() + 3 {
                    F::ZERO
                } else {
                    crate::lazy_constants::use_constants_table(
                        |table| metrics::mapped_c_sq(self, table),
                        table,
                        self.len() - 3,
                    )
                }
            }
        }
    }

    fn dc_sq_from_quadratic(&self) -> F {
        if self.len() != 4 {
            F::ZERO
        } else {
            let m_half = (-0.5_f32).into();
            let dv_0 = vector::sum_scaled(self, &[m_half, F::ONE, F::ZERO, m_half]);
            let dc2_0 = vector::length_sq(&dv_0);
            let dv_1 = vector::sum_scaled(self, &[m_half, F::ZERO, F::ONE, m_half]);
            let dc2_1 = vector::length_sq(&dv_1);
            utils::max(dc2_0, dc2_1)
        }
    }

    fn reduced_to_quadratic(&self) -> Option<Self> {
        if self.len() != 4 {
            None
        } else {
            Some(vec![
                self[0],
                vector::sum_scaled(
                    self,
                    &[
                        (-0.25_f32).into(),
                        0.75_f32.into(),
                        0.75_f32.into(),
                        (-0.25_f32).into(),
                    ],
                ),
                self[3],
            ])
        }
    }

    fn dc_sq_from_cubic(&self) -> F {
        todo!()
    }
    fn reduced_to_cubic(&self) -> Option<Self> {
        None
    }
}

impl<F: Num, const D: usize> BezierConstruct<F, D> for Vec<[F; D]> {
    fn of_builder(builder: &BezierBuilder<F, D>) -> Result<Self, ()> {
        let (mut matrix, pts) = builder.get_matrix_pts()?;
        let degree = pts.len() - 1;

        let mut lu = matrix.clone();
        let mut pivot = vec![0; degree + 1];
        let mut tr0 = vec![F::zero(); degree + 1];
        let mut tr1 = vec![F::zero(); degree + 1];

        if matrix::lup_decompose(degree + 1, &matrix, &mut lu, &mut pivot) == F::ZERO {
            return Err(());
        }

        assert!(
            matrix::lup_invert(degree + 1, &lu, &pivot, &mut matrix, &mut tr0, &mut tr1),
            "Matrix must be invertible"
        );
        let bezier: Vec<_> = matrix
            .chunks_exact(degree + 1)
            .map(|m| {
                m.iter()
                    .zip(pts.iter())
                    .fold([F::zero(); D], |acc, (coeff, p)| {
                        vector::add(acc, p, *coeff)
                    })
            })
            .collect();
        Ok(bezier)
    }
}
