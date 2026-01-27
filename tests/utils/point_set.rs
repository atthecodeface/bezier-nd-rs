use bezier_nd::BezierEval;
use bezier_nd::BezierSplit;

use bezier_nd::Float;
use bezier_nd::{Bezier, Num};
use geo_nd::{matrix, vector};

use super::{float_iter, max, min};

/// Create NPTS points on the bezier, and store them in a Vec
///
/// Split the bezier into lines given the straightness
///
/// For each line segment find NSEG points along the line and store all of these in a Vec
///
/// Find the largest closest distance between points on the bezier to segmented points
///
/// Find the largest closest distance between segmented points to points on the bezier
///
/// Both of these should be less than 'straightness'
pub struct BezierPtSet<F: Num, const D: usize> {
    pts: Vec<[F; D]>,
}

impl<F: Num, const D: usize> BezierPtSet<F, D> {
    pub fn is_empty(&self) -> bool {
        self.pts.is_empty()
    }
    pub fn len(&self) -> usize {
        self.pts.len()
    }
    pub fn of_point_at<B: BezierEval<F, [F; D]>>(bezier: &B, steps: usize) -> Self {
        let pts = float_iter(F::ZERO, F::ONE, steps)
            .map(|t| bezier.point_at(t))
            .collect();
        Self { pts }
    }

    pub fn of_point_iter<I: Iterator<Item = [F; D]>>(iter: I) -> Self {
        let pts = iter.collect();
        Self { pts }
    }

    pub fn of_line_iter<I: Iterator<Item = ([F; D], [F; D])>>(mut iter: I) -> Self {
        let l0 = iter.next().unwrap();
        let mut pts = vec![l0.0, l0.1];
        for (_, p) in iter {
            pts.push(p);
        }
        Self { pts }
    }

    pub fn of_line_iter_interpolated<I: Iterator<Item = ([F; D], [F; D])>>(
        iter: I,
        steps: usize,
    ) -> Self {
        let pts = iter
            .flat_map(|(p0, p1)| {
                float_iter(F::ZERO, F::ONE, steps).map(move |t| vector::mix(&p0, &p1, t))
            })
            .collect();
        Self { pts }
    }

    pub fn as_points(&self) -> impl Iterator<Item = [F; D]> + '_ {
        self.pts.iter().copied()
    }

    pub fn as_lines(&self) -> impl Iterator<Item = ([F; D], [F; D])> + '_ {
        self.pts
            .iter()
            .copied()
            .zip(self.pts.iter().skip(1).copied())
    }

    /// Calculate the maximum distance between successive points
    ///
    /// This helps to provide an error bound for a point supposedly 'on' the same Bezier
    pub fn max_pt_separation_sq(&self) -> F {
        let p0_p1_iter = self.pts.iter().zip(self.pts.iter().skip(1));
        p0_p1_iter.fold(F::ZERO, |max_d, (p0, p1)| {
            max(max_d, vector::distance_sq(p0, p1))
        })
    }

    #[track_caller]
    pub fn min_distance_sq_to_pt(&self, pt: &[F; D]) -> F {
        let mut min_d_sq = 1E10_f32.into();
        for p in &self.pts {
            let d_sq = vector::distance_sq(p, pt);
            min_d_sq = min(min_d_sq, d_sq);
        }
        min_d_sq
    }

    /// For BezierPtSet that are of the same length and corresponding points, find
    /// the maximum distance between any of the two corresponding points
    #[track_caller]
    pub fn max_distance_sq_between(&self, other: &Self) -> F {
        let mut max_d_sq = F::ZERO;
        assert_eq!(self.pts.len(), other.pts.len());
        for (p0, p1) in self.pts.iter().zip(other.pts.iter()) {
            let d_sq = vector::distance_sq(p0, p1);
            eprintln!("pts: {p0:?} {p1:?} d_sq:{d_sq}");
            if d_sq > max_d_sq {
                max_d_sq = d_sq;
            }
        }
        max_d_sq
    }

    /// For BezierPtSet for what should be (approx) the same Bezier's but with
    /// different t's for the points, find the max d_sq between *self*'s points and the closest of *other*'s points.
    ///
    /// This should be called only if other contains a full Bezier trace
    pub fn max_distance_sq_of_all_pts(&self, other: &Self) -> F {
        let mut max_excursion = F::zero();
        for p in self.pts.iter() {
            let min_d = other.pts.iter().fold(1.0E10_f32.into(), |md: F, q| {
                min(md, vector::distance_sq(p, q))
            });
            max_excursion = max(max_excursion, min_d);
        }
        max_excursion
    }
}
