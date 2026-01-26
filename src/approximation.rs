use crate::{BezierEval, BezierLineTIter, BezierPointTIter, BezierSplit, Num};
use geo_nd::vector;

pub struct Approximation<B, F, const D: usize>
where
    B: BezierEval<F, [F; D]> + BezierSplit + Clone,
    F: Num,
{
    bezier: B,
    straightness_sq: F,
    ts: Vec<F>,
    points: Vec<[F; D]>,
}

impl<B, F, const D: usize> Approximation<B, F, D>
where
    B: BezierEval<F, [F; D]> + BezierSplit + Clone,
    F: Num,
{
    pub fn new(bezier: &B, straightness_sq: F) -> Self {
        let points_and_ts: BezierPointTIter<_, _, _> =
            BezierLineTIter::new(bezier, straightness_sq).into();
        let mut points = vec![];
        let mut ts = vec![];
        for (t, p) in points_and_ts {
            points.push(p);
            ts.push(t);
        }
        let bezier = bezier.clone();
        Self {
            bezier,
            straightness_sq,
            ts,
            points,
        }
    }
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
            (find_it.unwrap_err(), false, false)
        }
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
                let t = (t - self.ts[n]) / dt;
                vector::sum_scaled(&self.points[n..n + 1], &[t, F::ONE - t])
            }
        }
    }
    fn derivative_at(&self, t: F) -> (F, [F; D]) {
        todo!()
    }
    fn endpoints(&self) -> (&[F; D], &[F; D]) {
        (&self.points.first().unwrap(), &self.points.last().unwrap())
    }
    fn closeness_sq_to_line(&self) -> F {
        F::ZERO
    }
    fn closeness_sq_to_quadratic(&self) -> F {
        F::ZERO
    }
    fn closeness_sq_to_cubic(&self) -> F {
        F::ZERO
    }
    fn num_control_points(&self) -> usize {
        self.points.len()
    }
    fn control_point(&self, n: usize) -> &[F; D] {
        &self.points[n]
    }
    fn degree(&self) -> usize {
        2
    }
}
