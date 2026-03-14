//a Imports
use crate::{BezierEval, BezierFlatIterator, BezierIterationType, Num};

/// This type provides optimized calculation of the distance from a point to a (set) of Bezier
///
/// The Beziers will be flattened to the given straightness, and if they form a closed path then
/// the winding will be used to determine if points are inside or outside the Bezier.
///
/// A clockwise winding has inside in the *right* of the lines that make up the curve; a
/// counter-clockwise winding has inside on the *left* of the lines that make up the curve.
// #[derive(Debug)]
pub struct FlattenedPath<F, P>
where
    F: Num,
    P: Clone,
    [P; 2]: BezierEval<F, P>,
{
    /// The straightness that the Bezier's will be flattened to
    closeness_sq: F,
    /// The winding, for a closed path, to determine inside/outside
    winding: bool,
    /// The beziers that make up the complete path
    beziers: Vec<Box<dyn BezierEval<F, P>>>,
    /// The points that form the flattened line segments
    flattened: Vec<(F, P)>,
    /// The indicies for each Bezier of the first (inclusive) and last (exclusive) flattened points
    bezier_indices: Vec<(usize, usize)>,
}

impl<F, P> FlattenedPath<F, P>
where
    F: Num,
    P: Clone,
    [P; 2]: BezierEval<F, P>,
{
    /// Create a new BezierDistance2D type given a straightness and winding
    ///
    /// For an open path, a winding of True suffices
    pub fn new(closeness_sq: F, winding: bool) -> Self {
        Self {
            closeness_sq,
            winding,
            beziers: vec![],
            flattened: vec![],
            bezier_indices: vec![],
        }
    }

    /// Add a bezier to the path
    pub fn add_bezier<B: BezierEval<F, P> + BezierFlatIterator<F, P> + Clone + 'static>(
        &mut self,
        bezier: &B,
    ) {
        self.beziers.push(Box::new(bezier.clone()));
        let first = self.flattened.len();
        for t_p in bezier.as_t_points(BezierIterationType::ClosenessSq(self.closeness_sq)) {
            self.flattened.push(t_p);
        }
        self.bezier_indices.push((first, self.flattened.len()));
    }

    /// Calculate the distance (squared) from the point to the closest point on the Bezier,
    /// and determine if the point is inside or outside the Beziers if they are a closed path
    pub fn distance_sq_to(&mut self, pt: &P) -> (F, bool) {
        let mut min_d2_s = (F::max_value(), 0);
        for (i, ((_t0, p0), (_t1, p1))) in self
            .flattened
            .iter()
            .zip(self.flattened.iter().skip(1))
            .enumerate()
        {
            let (_t, distance_sq) = [p0.clone(), p1.clone()].t_dsq_closest_to_pt(pt).unwrap();
            if distance_sq < min_d2_s.0 {
                min_d2_s = (distance_sq, i);
            }
        }
        (min_d2_s.0, false)
    }
}
