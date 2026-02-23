use crate::{
    bernstein_fns, constants, metrics, BezierBuilder, BezierConstruct, BezierElevate, BezierError,
    BezierEval, BezierFlatIterator, BezierIterationType, BezierLineTIter, BezierMap, BezierMetric,
    BezierOps, BezierReduce, BezierReduction, BezierSplit, Num, CONSTANTS_TABLE,
};

use geo_nd::matrix;
use geo_nd::vector;

/// This type stores a Bernstein Bezier of up to `N` control points each of
/// dimension `D`, with coordinate/parameter values of type `F`
///
/// This is deliberately [Cpoy], and contains an array of points. Only the first
/// `degree+1` points are used. As this supports arbitrary degree Bezier curves
/// (up to degree `N-1`), it can provide trait implementations of
/// [BezierReduce], [BezierElevate], and [BezierMapped] which produce `Self` as
/// the resultant type.
///
/// If truly arbitrary degree Bezier curves are required then a `Vec<F>` can be
/// used, but this obviously requires `Alloc`.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct BezierND<F, const N: usize, const D: usize>
where
    F: Num,
{
    /// Degree, which is one less than number of valid control points
    ///
    /// i.e. for a cubic this is three
    ///
    /// A Bezier of degree 0 is a fixed point
    degree: usize,
    /// Control points - 0..=degree are valid
    pts: [[F; D]; N],
}

impl<F, const N: usize, const D: usize> std::default::Default for BezierND<F, N, D>
where
    F: Num,
{
    fn default() -> Self {
        let pts = [[F::ZERO; D]; N];
        Self { degree: 0, pts }
    }
}

impl<F, const N: usize, const D: usize> std::fmt::Display for BezierND<F, N, D>
where
    F: Num,
{
    //mp fmt - format a `Bezier` for display
    /// Display the `Bezier' as sets of points
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[")?;
        let mut needs_comma = false;
        for p in self.pts.iter().take(self.degree + 1) {
            if needs_comma {
                write!(f, ", ")?;
            }
            needs_comma = true;
            vector::fmt(f, p)?;
        }
        write!(f, "]")
    }

    //zz All done
}

//ip Bezier iterators
impl<F, const N: usize, const D: usize> BezierND<F, N, D>
where
    F: Num,
{
    /// Find the maximum degree this type can encode
    pub fn max_degree() -> usize {
        N - 1
    }

    //mp Reduce-and-split iterator
    /// Apply a (new_degree+1) by (degree+1) matrix to the points to generate a new Bezier
    /// of a new degree
    pub fn reduce_and_split_iter<'a>(
        &'a self,
        reduce_matrix: &'a [F],
        elev_reduce_matrix: &'a [F],
        reduce_degree: usize,
        max_dc_sq: F,
    ) -> BezierReduceIter<'a, F, N, D> {
        BezierReduceIter {
            reduce_matrix,
            elev_reduce_matrix,
            reduce_degree,
            max_dc_sq,
            stack: vec![],
        }
    }
}

//tp BezierReduceIter
/// A type that provided an iterator implementaion for splitting a Bezier
/// into reduced Beziers given a maximum 'dc' metric
#[derive(Clone, PartialEq, Debug)]
pub struct BezierReduceIter<'a, F, const N: usize, const D: usize>
where
    F: Num,
{
    reduce_matrix: &'a [F],
    elev_reduce_matrix: &'a [F],
    reduce_degree: usize,
    max_dc_sq: F,
    stack: Vec<(usize, BezierND<F, N, D>)>,
}

impl<F, const N: usize, const D: usize> BezierND<F, N, D>
where
    F: Num,
{
    /// Create a new Bezier given all the control points
    pub fn new(pts: &[[F; D]]) -> Self {
        let n = pts.len();
        assert!(n >= 1, "Beziers have at least 1 point");
        assert!(
            n <= N,
            "Attempt to create a Bezier of max {N} pts with actually {n} pts"
        );
        let mut s = Self {
            degree: n - 1,
            ..Default::default()
        };
        s.pts.split_at_mut(n).0.copy_from_slice(pts);
        s
    }

    /// Create a new Bezier that is the nth derivative of this
    /// subject to a scaling factor (i.e. the actual Bezier should be scaled up by F)
    ///
    /// As this type uses Bernstein polynomials, this is the Bezier of degree N-n
    /// that has points control points:
    ///
    /// `Pn[i] = (N n) . Sum((-1)^j * P[i+n-j] * (n j))`
    ///
    /// (n j) is kept in `bezier_nd::constants::BINOMIALS[n][j+1]`
    pub fn nth_derivative(&self, n: usize) -> (Self, F) {
        let mut s = Self {
            degree: self.degree - n,
            ..Default::default()
        };
        s.degree = self.degree - n;
        let scale = bernstein_fns::nth_derivative(&self.pts[0..self.degree + 1], n, &mut s.pts);
        (s, scale)
    }

    /// Apply a (new_degree+1) by (degree+1) matrix to the points to generate a new Bezier
    /// of a new degree
    pub fn apply_matrix(&self, matrix: &[F], new_degree: usize) -> Self {
        assert!(
            new_degree < N + 1,
            "Cannot create a Bezier<{N}> of {new_degree}",
        );
        assert_eq!(
            matrix.len(),
            (new_degree + 1) * (self.degree + 1),
            "Matrix to apply to Bezier of degree {} to degree {new_degree} must have {} elements",
            self.degree,
            (new_degree + 1) * (self.degree + 1)
        );
        let _nr = new_degree + 1;
        let nc = self.degree + 1;
        let mut s = Self {
            degree: new_degree,
            ..Default::default()
        };
        for (m, sp) in matrix.chunks_exact(nc).zip(s.pts.iter_mut()) {
            let mut sum = [F::zero(); D];
            for (coeff, p) in m.iter().zip(self.pts.iter()) {
                sum = vector::add(sum, p, *coeff);
            }
            *sp = sum;
        }
        s
    }

    /// Return a slice of the control points
    pub fn pts(&self) -> &[[F; D]] {
        &self.pts[0..self.degree + 1]
    }

    /// Return a mutable slice of the control points
    pub fn pts_mut(&mut self) -> &mut [[F; D]] {
        &mut self.pts[0..self.degree + 1]
    }

    /// Apply a (degree+1) by (degree+1) matrix (should be elevate-of-reduce) to the points
    /// and calculate the new dc squared
    pub fn dc2_of_ele_red(&self, matrix: &[F]) -> F {
        assert_eq!(
            matrix.len(),
            (self.degree + 1) * (self.degree + 1),
            "Matrix to apply to Bezier of degree {} must have {} elements",
            self.degree,
            (self.degree + 1) * (self.degree + 1)
        );
        let mut max_d2 = F::zero();
        let nc = self.degree + 1;
        for (m, s) in matrix.chunks_exact(nc).zip(self.pts.iter()) {
            let mut sum = [F::zero(); D];
            for (coeff, p) in m.iter().zip(self.pts.iter()) {
                sum = vector::add(sum, p, *coeff);
            }
            let d2 = vector::distance_sq(s, &sum);
            if max_d2 < d2 {
                max_d2 = d2;
            }
        }
        max_d2
    }
}

impl<F, const N: usize, const D: usize> BezierEval<F, [F; D]> for BezierND<F, N, D>
where
    F: Num,
{
    fn distance_sq_between(&self, p0: &[F; D], p1: &[F; D]) -> F {
        vector::distance_sq(p0, p1)
    }
    fn point_at(&self, t: F) -> [F; D] {
        bernstein_fns::values::point_at(self.pts(), t)
    }
    fn derivative_at(&self, t: F) -> (F, [F; D]) {
        bernstein_fns::values::derivative_at(self.pts(), t)
    }

    fn endpoints(&self) -> ([F; D], [F; D]) {
        (self.pts[0], self.pts[self.degree])
    }

    fn closeness_sq_to_line(&self) -> F {
        metrics::dc_sq_from_line(&self.pts[0..=self.degree])
    }

    fn dc_sq_from_line(&self) -> F {
        metrics::dc_sq_from_line(&self.pts[0..=self.degree])
    }

    fn num_control_points(&self) -> usize {
        self.degree + 1
    }

    fn control_points(&self) -> &[[F; D]] {
        &self.pts[0..=self.degree]
    }

    fn degree(&self) -> usize {
        self.degree
    }
    fn metric_from(&self, other: Option<&[[F; D]]>, metric: BezierMetric) -> Option<F> {
        if let Some(other) = other {
            metrics::metric_from(&self.pts[0..=self.degree], other, metric)
        } else {
            Some(metrics::metric_from_line(
                &self.pts[0..=self.degree],
                metric,
            ))
        }
    }

    fn t_coords_at_min_max(
        &self,
        pt_index: usize,
        give_min: bool,
        give_max: bool,
    ) -> (Option<(F, F)>, Option<(F, F)>) {
        // Generally defer to the methods for Beziers of the associated degree
        //
        // When the value is not analytically calculable return None
        match self.degree {
            0 => (
                give_min.then_some((F::ZERO, self.pts[0][pt_index])),
                give_max.then_some((F::ZERO, self.pts[0][pt_index])),
            ),
            1 => <&[[F; D]] as TryInto<&[[F; D]; 2]>>::try_into(&self.pts[0..2])
                .unwrap()
                .t_coords_at_min_max(pt_index, give_min, give_max),
            2 => <&[[F; D]] as TryInto<&[[F; D]; 3]>>::try_into(&self.pts[0..3])
                .unwrap()
                .t_coords_at_min_max(pt_index, give_min, give_max),

            3 => <&[[F; D]] as TryInto<&[[F; D]; 4]>>::try_into(&self.pts[0..4])
                .unwrap()
                .t_coords_at_min_max(pt_index, give_min, give_max),

            _ => (None, None),
        }
    }

    fn t_dsq_closest_to_pt(&self, pt: &[F; D]) -> Option<(F, F)> {
        metrics::t_dsq_closest_to_pt(self.pts(), pt)
    }

    fn est_min_distance_sq_to(&self, pt: &[F; D]) -> F {
        if self.pts.len() == 0 {
            F::ZERO
        } else {
            metrics::est_min_distance_sq_to(self.pts(), pt)
        }
    }
}

impl<F, const N: usize, const D: usize> BezierOps<F, [F; D]> for BezierND<F, N, D>
where
    F: Num,
{
    fn add(&mut self, other: &Self) -> bool {
        if self.degree != other.degree {
            false
        } else {
            for (s, o) in self.pts_mut().iter_mut().zip(other.pts.iter()) {
                *s = vector::add(*s, o, F::ONE)
            }
            true
        }
    }

    /// Subtract another Bezier from this
    fn sub(&mut self, other: &Self) -> bool {
        if self.degree != other.degree {
            false
        } else {
            for (s, o) in self.pts_mut().iter_mut().zip(other.pts.iter()) {
                *s = vector::sub(*s, o, F::ONE)
            }
            true
        }
    }

    /// Scale the points
    fn scale(&mut self, scale: F) {
        self.map_pts(&|_, p| vector::scale(*p, scale));
    }

    /// Map the points
    fn map_pts(&mut self, map: &dyn Fn(usize, &[F; D]) -> [F; D]) {
        for (i, p) in self.pts.iter_mut().take(self.degree + 1).enumerate() {
            *p = map(i, p);
        }
    }
    fn map_all_pts<'a>(&'a mut self, map: &'a mut dyn FnMut(&'a mut [[F; D]]) -> bool) -> bool {
        map(&mut self.pts[0..self.degree + 1])
    }
}

impl<F, const N: usize, const D: usize> BezierSplit<F> for BezierND<F, N, D>
where
    F: Num,
{
    fn split(&self) -> (Self, Self) {
        self.split_at(F::frac(1, 2))
    }
    fn split_at(&self, t: F) -> (Self, Self) {
        let mut first = *self;
        let mut latter = *self;
        bernstein_fns::de_casteljau::split_at(latter.pts_mut(), t, first.pts_mut());
        (first, latter)
    }
    fn section(&self, t0: F, t1: F) -> Self {
        let mut to_split = *self;
        if t0 > F::ZERO {
            bernstein_fns::de_casteljau::bezier_from(to_split.pts_mut(), t0);
        }
        if t1 < F::ONE {
            let t10 = (t1 - t0) / (F::ONE - t0);
            bernstein_fns::de_casteljau::bezier_to(to_split.pts_mut(), t10);
        }
        to_split
    }
}

impl<F, const N: usize, const D: usize> BezierFlatIterator<F, [F; D]> for BezierND<F, N, D>
where
    F: Num,
{
    fn as_t_lines(
        &self,
        iter_type: BezierIterationType<F>,
    ) -> impl Iterator<Item = (F, [F; D], F, [F; D])> {
        BezierLineTIter::of_iter_type(self, iter_type)
    }
}

impl<F: Num, const N: usize, const D: usize> BezierConstruct<F, D> for BezierND<F, N, D> {
    fn of_degree(degree: usize) -> Result<Self, BezierError> {
        if degree + 1 < N {
            let mut s = Self::default();
            s.degree = degree;
            Ok(s)
        } else {
            Err(BezierError::BadBuildDegree(degree))
        }
    }

    fn of_builder(builder: &BezierBuilder<F, D>) -> Result<Self, BezierError> {
        let (mut matrix, pts) = builder.get_matrix_pts()?;
        if pts.len() > N {
            return Err(BezierError::MaxBuildDegree(N - 1, pts.len() + 1));
        }
        let degree = pts.len() - 1;
        let bezier = Self::new(&pts);

        let mut lu = matrix.clone();
        let mut pivot = vec![0; degree + 1];
        let mut tr0 = vec![F::zero(); degree + 1];
        let mut tr1 = vec![F::zero(); degree + 1];

        if matrix::lup_decompose(degree + 1, &matrix, &mut lu, &mut pivot) == F::ZERO {
            return Err(BezierError::BadBuildConstraints);
        }

        assert!(
            matrix::lup_invert(degree + 1, &lu, &pivot, &mut matrix, &mut tr0, &mut tr1),
            "Matrix must be invertible"
        );
        Ok(bezier.apply_matrix(&matrix, degree))
    }
}

//ip Iterator for BezierReduceIter
impl<'a, F, const N: usize, const D: usize> std::iter::Iterator for BezierReduceIter<'a, F, N, D>
where
    F: Num,
{
    type Item = (usize, BezierND<F, N, D>);
    fn next(&mut self) -> Option<(usize, BezierND<F, N, D>)> {
        while let Some((n, b)) = self.stack.pop() {
            let dc2 = b.dc2_of_ele_red(self.elev_reduce_matrix);
            if dc2 > self.max_dc_sq {
                let (b0, b1) = b.split();
                self.stack.push((n + 1, b1));
                self.stack.push((n + 1, b0));
            } else {
                return Some((n, b.apply_matrix(self.reduce_matrix, self.reduce_degree)));
            }
        }
        None
    }
}
impl<F, const N: usize, const D: usize> BezierElevate<F, [F; D]> for BezierND<F, N, D>
where
    F: crate::Num,
{
    type ElevatedByOne = Self;

    fn elevate_by_one(&self) -> Option<Self> {
        if self.degree + 1 >= N {
            None
        } else {
            // n = number of points, i.e. self.pts[n-1] is the last valid point
            let n = self.degree + 1;
            let mut s = Self {
                degree: n,
                ..Default::default()
            };
            s.pts[0] = self.pts[0];
            s.pts[s.degree] = self.pts[self.degree];
            let n_f = F::of_usize(n);
            let mut s0 = F::ONE;
            for (p, (p0, p1)) in s
                .pts
                .iter_mut()
                .skip(1)
                .zip(self.pts.iter().zip(self.pts.iter().skip(1)))
            {
                *p = vector::sum_scaled(&[*p0, *p1], &[s0 / n_f, F::ONE - s0 / n_f]);
                s0 += F::ONE;
            }
            Some(s)
        }
    }
}

impl<F, const N: usize, const D: usize> BezierReduce<F, [F; D]> for BezierND<F, N, D>
where
    F: crate::Num,
{
    type Reduced = Self;

    fn can_reduce(&self, method: BezierReduction) -> bool {
        constants::reduce_table_of_match(method, self.degree).is_some()
    }
    fn reduce(&self, method: BezierReduction) -> Option<Self::Reduced> {
        if let Some(table) = constants::reduce_table_of_match(method, self.degree) {
            let mut result = Self::default();
            CONSTANTS_TABLE.use_constants_table(
                |table| {
                    bernstein_fns::transform::transform_pts(
                        table,
                        &self.pts[0..=self.degree],
                        &mut result.pts[0..self.degree],
                    )
                },
                table,
            );
            result.degree = self.degree - 1;
            Some(result)
        } else {
            None
        }
    }
    fn dc_sq_from_reduction(&self, method: BezierReduction) -> F {
        if let Some(table) = constants::er_minus_i_table_of_match(method, self.degree) {
            CONSTANTS_TABLE.use_constants_table(
                |table| metrics::mapped_c_sq(&self.pts[0..=self.degree], table),
                table,
            )
        } else {
            F::ZERO
        }
    }
}

impl<F, const N: usize, const D: usize> BezierMap<F, [F; D]> for BezierND<F, N, D>
where
    F: crate::Num,
{
    type Mapped = Self;

    fn mapped_to_degree(&self, to_degree: usize, matrix: &[F]) -> Option<Self::Mapped> {
        if to_degree + 1 >= N {
            return None;
        }
        assert_eq!(
            matrix.len(),
            (to_degree + 1) * (self.degree+1),
            "Matrix for mapping from degree {} to degree {to_degree} must be {}x{} but was {} in total",
            self.degree,
            to_degree + 1,
            self.degree+1,
            matrix.len()
        );
        let mut result = Self::default();
        result.degree = to_degree;
        bernstein_fns::transform::transform_pts(matrix, &self.pts, &mut result.pts);
        Some(result)
    }

    fn dc_sq_of_mapped_from_line(&self, to_degree: usize, matrix: &[F]) -> Option<F> {
        if to_degree + 1 >= N {
            return None;
        }
        assert_eq!(
            matrix.len(),
            (to_degree + 1) * (self.degree+1),
            "Matrix for mapping from degree {} to degree {to_degree} must be {}x{} but was {} in total",
            self.degree,
            to_degree + 1,
            self.degree+1,
            matrix.len()
        );
        Some(metrics::dc_sq_mapped_from_line(
            &self.pts, to_degree, matrix,
        ))
    }
}
