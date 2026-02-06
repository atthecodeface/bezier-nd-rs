use crate::{BezierEval, BezierIntoIterator, BezierSplit, Float, Num};
use geo_nd::vector;

/// An approximation to a Bezier of potentially high degree,
/// consisting of an array of points that join to form line
/// segments that approximate the Bezier within a closeness_sq
/// such that all points on the line segments are within that closeness_sq
/// of the corresponding Bezier point
///
/// The type needs to contain a clone of the original Bezier,
/// to provide correct values for some of the evaulation methods.
///
/// This type should provide implementations of
///
/// BezierEval
/// BezierDistance
/// BezierMinMax
#[derive(Debug, Clone)]
pub struct Approximation<B, F, const D: usize>
where
    B: BezierEval<F, [F; D]> + BezierSplit + Clone,
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
    B: BezierEval<F, [F; D]> + BezierSplit + Clone,
    F: Num,
{
    /// Create a new Approximation from a Bezier and a given closeness_sq
    ///
    /// This maintains the mapping of `t` on the original Bezier to the
    /// approximation, and so may use more line segments than are required
    /// just to approximate for simple straightness
    pub fn new(bezier: &B, closeness_sq: F) -> Self {
        let points_and_ts = bezier.as_t_points_dc(closeness_sq);
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

impl<B, F, const D: usize> Approximation<B, F, D>
where
    B: BezierEval<F, [F; D]> + BezierSplit + Clone,
    F: Float,
{
    /// Calculates the length of a setion of the Bezier
    pub fn section_length(&self, mut t0: F, mut t1: F) -> F {
        if t0 < F::ZERO {
            t0 = F::ZERO;
        }
        if t1 > F::ONE {
            t1 = F::ONE;
        }

        if t0 > t1 {
            self.section_length(t1, t0)
        } else if t0 == t1 {
            F::ZERO
        } else {
            let (pt_0, _, _) = self.find_pts(t0);
            let (pt_1, _, _) = self.find_pts(t1);
            let mut t = t0;
            let mut total = F::ZERO;
            for i in pt_0..=pt_1 {
                if i + 1 >= self.points.len() {
                    break;
                }
                if self.ts[i] > t1 {
                    break;
                }
                if t < self.ts[i] {
                    t = self.ts[i];
                }
                let dt = self.ts[i + 1] - self.ts[i];
                let fract = (self.ts[i + 1] - t) / dt;
                total += vector::distance(&self.points[i], &self.points[i + 1]) * fract
            }
            total
        }
    }

    /// Calculates the length of a setion of the Bezier
    pub fn t_of_distance(&self, mut distance: F) -> Option<F> {
        for (i, (p0, p1)) in self
            .points
            .iter()
            .zip(self.points.iter().skip(1))
            .enumerate()
        {
            let dp = vector::distance(p0, p1);
            if distance <= dp {
                return Some(self.ts[i] + (self.ts[1] - self.ts[0]) * distance / dp);
            }
            distance -= dp;
        }
        None
    }
}

impl<B, F, const D: usize> BezierEval<F, [F; D]> for Approximation<B, F, D>
where
    B: BezierEval<F, [F; D]> + BezierSplit + Clone,
    F: Num,
{
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
                eprintln!("{t} : {t_rel}, {n} {t0} {t1} {dt} {t_rel}");
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
    fn endpoints(&self) -> (&[F; D], &[F; D]) {
        (self.points.first().unwrap(), self.points.last().unwrap())
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
    fn control_point(&self, n: usize) -> &[F; D] {
        self.bezier.control_point(n)
    }
    fn degree(&self) -> usize {
        self.bezier.degree()
    }
    fn for_each_control_points(&self, map: &mut dyn FnMut(&[F; D])) {
        self.bezier.for_each_control_points(map);
    }
}

impl<B, F, const D: usize> BezierIntoIterator<F, [F; D]> for Approximation<B, F, D>
where
    B: BezierEval<F, [F; D]> + BezierSplit + Clone,
    F: Num,
{
    fn as_lines(&self, _closeness_sq: F) -> impl Iterator<Item = ([F; D], [F; D])> {
        self.points
            .iter()
            .copied()
            .zip(self.points.iter().copied().skip(1))
    }

    fn as_points(&self, _closeness_sq: F) -> impl Iterator<Item = [F; D]> {
        self.points.iter().copied()
    }

    fn as_t_lines(&self, _closeness_sq: F) -> impl Iterator<Item = (F, [F; D], F, [F; D])> {
        self.iter_t_pts()
            .zip(self.iter_t_pts().skip(1))
            .map(|((a, b), (c, d))| (a, b, c, d))
    }

    fn as_t_points(&self, _closeness_sq: F) -> impl Iterator<Item = (F, [F; D])> {
        self.iter_t_pts()
    }

    fn as_t_lines_dc(&self, closeness_sq: F) -> impl Iterator<Item = (F, [F; D], F, [F; D])> {
        self.as_t_lines(closeness_sq)
    }
    fn as_t_points_dc(&self, _closeness_sq: F) -> impl Iterator<Item = (F, [F; D])> {
        self.iter_t_pts()
    }
}
