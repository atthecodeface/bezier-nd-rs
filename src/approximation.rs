use crate::{
    utils, BezierEval, BezierFlatIterator, BezierIterationType, BezierMetric, BezierSplit, Num,
};
use geo_nd::vector;

/// An approximation to a Bezier of potentially high degree,
/// consisting of an array of points that join to form line
/// segments that approximate the Bezier within a closeness_sq
/// such that all points on the line segments are within that closeness_sq
/// of the corresponding Bezier point
///
/// The type needs to contain a clone of the original Bezier,
/// to provide correct values for some of the evaulation methods.
#[derive(Debug, Clone)]
pub struct Approximation<B, F, const D: usize>
where
    B: BezierEval<F, [F; D]> + BezierSplit<F> + BezierFlatIterator<F, [F; D]> + Clone,
    F: Num,
{
    /// The underlying Bezier
    bezier: B,

    /// The closeness_sq for each point on the approximation to
    /// the associated point on the Bezier
    closeness_sq: F,

    /// The array of 't' parameter values for the original Bezier that
    /// the approximation has for its points.
    ///
    /// The first of these will always be 0, and the last 1
    ts: Vec<F>,

    /// The array of points (corresponding to the 't' parameter values)
    /// that make up the approximation
    ///
    /// The first of these will be the starting point of the original Bezier,
    /// and the last will be the ending point of the original Bezier
    points: Vec<[F; D]>,
}

impl<B, F, const D: usize> Approximation<B, F, D>
where
    B: BezierEval<F, [F; D]> + BezierSplit<F> + BezierFlatIterator<F, [F; D]> + Clone,
    F: Num,
{
    /// Create a new Approximation from a Bezier and a given closeness_sq
    ///
    /// This maintains the mapping of `t` on the original Bezier to the
    /// approximation, and so may use more line segments than are required
    /// just to approximate for simple straightness
    pub fn new(bezier: &B, closeness_sq: F) -> Self {
        let points_and_ts = bezier.as_t_points(BezierIterationType::DcClosenessSq(closeness_sq));
        let mut points = vec![];
        let mut ts = vec![];
        for (t, p) in points_and_ts {
            points.push(p);
            ts.push(t);
        }
        let bezier = bezier.clone();
        Self {
            bezier,
            closeness_sq,
            ts,
            points,
        }
    }

    /// Create a new Approximation from a Bezier and regular `t` intervals
    ///
    /// The specified `closeness_sq` should *exceed* the expected distance; this
    /// can be taken (if required) from metrics::dc_sq_from_line(bezier)
    pub fn of_regular_t(bezier: &B, num_pts: usize, closeness_sq: Option<F>) -> Self {
        let ts: Vec<F> = utils::float_iter(num_pts).collect();
        let points = ts.iter().map(|t| bezier.point_at(*t)).collect();
        let closeness_sq = closeness_sq.unwrap_or_else(|| bezier.dc_sq_from_line());
        let bezier = bezier.clone();
        Self {
            bezier,
            closeness_sq,
            ts,
            points,
        }
    }

    /// Find the index into the Approximation of the point/t value that
    /// is at or just before the given t value; indicate also if this matches
    /// a t value precisely.
    ///
    /// The return value for a t value that is not precisely in the approximation
    /// is (index, false, false) where index has a t value is less than that requested,
    /// and where the next index has a t value that is greater than that requested.
    ///
    /// The return value for a t value that is precisely the same as a point in the
    /// approximation, but t<1, is (index, true, false), where index has the matching t value
    ///
    /// The return value for t=1.0 is (index, false, true), where index is the last point
    /// (and hence must have t=1.0)
    pub fn find_pts(&self, t: F) -> (usize, bool, bool) {
        assert!(t >= F::ZERO);
        assert!(t <= F::ONE);
        let find_it = self.ts.binary_search_by(|pt| pt.partial_cmp(&t).unwrap());
        if let Ok(n) = find_it {
            if t < F::ONE {
                (n, true, false)
            } else {
                (n, false, true)
            }
        } else {
            (find_it.unwrap_err() - 1, false, false)
        }
    }

    /// Iterate over the t values and points in the approximation
    ///
    /// This is similar to the 'as_t_points' method on Beziers (when implemented)
    /// with the 'closeness_sq' having been used in the generation of the Approximation
    pub fn iter_t_pts(&self) -> impl Iterator<Item = (F, [F; D])> + '_ {
        self.ts
            .iter()
            .zip(self.points.iter())
            .map(|(a, b)| (*a, *b))
    }

    /// Borrow the points of the approximation
    pub fn points(&self) -> &[[F; D]] {
        &self.points
    }

    /// Borrow the t values of the approximation
    pub fn ts(&self) -> &[F] {
        &self.ts
    }

    /// Iterate over the lines
    pub fn iter_lines(&self) -> impl Iterator<Item = (&[F; D], &[F; D])> + '_ {
        self.points.iter().zip(self.points.iter().skip(1))
    }
}

impl<B, F, const D: usize> BezierEval<F, [F; D]> for Approximation<B, F, D>
where
    B: BezierEval<F, [F; D]> + BezierSplit<F> + BezierFlatIterator<F, [F; D]> + Clone,
    F: Num,
{
    fn distance_sq_between(&self, p0: &[F; D], p1: &[F; D]) -> F {
        vector::distance_sq(p0, p1)
    }
    fn point_at(&self, t: F) -> [F; D] {
        if t <= F::ZERO {
            self.points[0]
        } else if t >= F::ONE {
            *self.points.last().unwrap()
        } else {
            let (n, is_start, is_end) = self.find_pts(t);
            if is_start || is_end {
                self.points[n]
            } else {
                let t0 = self.ts[n];
                let t1 = self.ts[n + 1];
                let dt = t1 - t0;
                let t_rel = (t - t0) / dt;
                vector::sum_scaled(&self.points[n..], &[F::ONE - t_rel, t_rel])
            }
        }
    }
    fn derivative_at(&self, t: F) -> (F, [F; D]) {
        if t < F::ZERO || t > F::ONE {
            (F::ZERO, self.points[0])
        } else {
            let (n, _is_start, is_end) = self.find_pts(t);
            let n = if is_end { n - 1 } else { n };
            let t0 = self.ts[n];
            let t1 = self.ts[n + 1];
            let dt = t1 - t0;
            let dt_r = F::ONE / dt;
            (
                F::ONE,
                vector::sum_scaled(&self.points[n..], &[-dt_r, dt_r]),
            )
        }
    }
    fn endpoints(&self) -> ([F; D], [F; D]) {
        (*self.points.first().unwrap(), *self.points.last().unwrap())
    }
    fn closeness_sq_to_line(&self) -> F {
        self.bezier.closeness_sq_to_line() + self.closeness_sq
    }
    fn dc_sq_from_line(&self) -> F {
        self.bezier.dc_sq_from_line() + self.closeness_sq
    }
    fn num_control_points(&self) -> usize {
        self.bezier.num_control_points()
    }
    fn control_points(&self) -> &[[F; D]] {
        self.bezier.control_points()
    }
    fn degree(&self) -> usize {
        self.bezier.degree()
    }
    fn metric_from(&self, other: Option<&[[F; D]]>, metric: BezierMetric) -> Option<F> {
        self.bezier.metric_from(other, metric)
    }
    // Can probably do better
    fn t_dsq_closest_to_pt(&self, pt: &[F; D]) -> Option<(F, F)> {
        self.bezier.t_dsq_closest_to_pt(pt)
    }
    fn est_min_distance_sq_to(&self, _p: &[F; D]) -> F {
        F::ZERO
    }
    fn t_coords_at_min_max(
        &self,
        pt_index: usize,
        give_min: bool,
        give_max: bool,
    ) -> (Option<(F, F)>, Option<(F, F)>) {
        assert!(
            pt_index < D,
            "point index out of range of point array size {D}"
        );
        let mut min: Option<(F, F)> = None;
        let mut max: Option<(F, F)> = None;
        for (t, p) in self.iter_t_pts() {
            let p = p[pt_index];
            if give_min {
                if let Some(m) = min {
                    if p < m.1 {
                        min = Some((t, p));
                    }
                } else {
                    min = Some((t, p));
                }
            }
            if give_max {
                if let Some(m) = max {
                    if p > m.1 {
                        max = Some((t, p));
                    }
                } else {
                    max = Some((t, p));
                }
            }
        }
        (min, max)
    }
}

impl<B, F, const D: usize> BezierFlatIterator<F, [F; D]> for Approximation<B, F, D>
where
    B: BezierEval<F, [F; D]> + BezierSplit<F> + BezierFlatIterator<F, [F; D]> + Clone,
    F: Num,
{
    fn as_t_lines(
        &self,
        _iter_type: BezierIterationType<F>,
    ) -> impl Iterator<Item = (F, [F; D], F, [F; D])> {
        self.iter_t_pts()
            .zip(self.iter_t_pts().skip(1))
            .map(|((a, b), (c, d))| (a, b, c, d))
    }
}
