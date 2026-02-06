use crate::Num;
use crate::{
    bernstein_fns, BezierBuilder, BezierConstruct, BezierElevate, BezierEval, BezierMinMax,
    BezierOps, BezierSection, BezierSplit,
};
use geo_nd::{matrix, vector};

impl<F: Num, const D: usize> BezierEval<F, [F; D]> for Vec<[F; D]> {
    fn point_at(&self, t: F) -> [F; D] {
        let degree = self.len() - 1;

        self.iter()
            .zip(bernstein_fns::basis_coeff_enum_num(degree, t))
            .fold([F::ZERO; D], |acc, (pt, (_i, coeff))| {
                vector::add(acc, pt, coeff)
            })
    }
    fn derivative_at(&self, t: F) -> (F, [F; D]) {
        let degree = self.len() - 1;
        let (reduce, coeffs) = bernstein_fns::basis_dt_coeff_enum_num(degree, t);
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
    fn for_each_control_points(&self, map: &mut dyn FnMut(&[F; D])) {
        self.iter().for_each(map)
    }
}

impl<F: Num, const D: usize> BezierOps<F, [F; D]> for Vec<[F; D]> {
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

impl<F: Num, const D: usize> BezierMinMax<F> for Vec<[F; D]> {
    fn t_coord_at_min_max(&self, use_max: bool, pt_index: usize) -> Option<(F, F)> {
        match self.len() {
            1 => Some((F::ZERO, self[0][pt_index])),
            2 => <&[[F; D]] as TryInto<&[[F; D]; 2]>>::try_into(&self[0..2])
                .unwrap()
                .t_coord_at_min_max(use_max, pt_index),
            3 => <&[[F; D]] as TryInto<&[[F; D]; 3]>>::try_into(&self[0..3])
                .unwrap()
                .t_coord_at_min_max(use_max, pt_index),
            4 => <&[[F; D]] as TryInto<&[[F; D]; 4]>>::try_into(&self[0..4])
                .unwrap()
                .t_coord_at_min_max(use_max, pt_index),
            _ => None,
        }
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

impl<F: 'static + Num, const D: usize> BezierElevate<F, [F; D]> for Vec<[F; D]> {
    type ElevatedByOne = Self;
    type Elevated = Self;
    fn elevate_by_one(&self) -> Option<Self> {
        let mut s = Vec::with_capacity(self.len() + 1);
        s.push(*self.first().unwrap());
        let n = self.len();
        let n_f: F = (n as f32).into();
        for (i, (p0, p1)) in self.iter().zip(self.iter().skip(1)).enumerate() {
            let s0: F = ((i + 1) as f32).into();
            let s1: F = ((n - (i + 1)) as f32).into();
            s.push(vector::sum_scaled(&[*p0, *p1], &[s0 / n_f, s1 / n_f]));
        }
        s.push(*self.last().unwrap());
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
