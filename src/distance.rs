use std::cmp::min;

//a Imports
use geo_nd::vector;
use geo_nd::Float;

use crate::{Bezier, BezierLineIter, BezierPointIter};

/// This type provides optimized calculation of the distance from a point to a (set) of Bezier
///
/// The Beziers will be flattened to the given straightness, and if they form a closed path then
/// the winding will be used to determine if points are inside or outside the Bezier.
///
/// A clockwise winding has inside in the *right* of the lines that make up the curve; a
/// counter-clockwise winding has inside on the *left* of the lines that make up the curve.
pub struct BezierDistance2D<F: Float> {
    /// The straightness that the Bezier's will be flattened to
    straightness: F,
    /// The winding, for a closed path, to determine inside/outside
    winding: bool,
    /// The beziers that make up the complete path
    beziers: Vec<Bezier<F, 2>>,
    /// The flattened line segments
    flattened: Vec<[F; 2]>,
    /// The indicies for each Bezier of the first and last flattened line segments
    bezier_indices: Vec<(usize, usize)>,
}

impl<F: Float> BezierDistance2D<F> {
    /// Create a new BezierDistance2D type given a straightness and winding
    ///
    /// For an open path, a winding of True suffices
    pub fn new(straightness: F, winding: bool) -> Self {
        Self {
            straightness,
            winding,
            beziers: vec![],
            flattened: vec![],
            bezier_indices: vec![],
        }
    }

    /// Add a bezier to the path
    pub fn add_bezier(&mut self, bezier: Bezier<F, 2>) {
        assert!(
            self.flattened.is_empty(),
            "Addind Beziers is only permitted to a new or reset BezierDistance2D"
        );
        self.beziers.push(bezier.clone());
    }

    /// Add a bezier to the path, splitting it into smaller beziers to get within a straightness
    pub fn add_bezier_with_straightness(&mut self, bezier: Bezier<F, 2>, straightness: F) {
        if bezier.is_straight(straightness) {
            self.beziers.push(bezier);
        } else {
            let (b0, b1) = bezier.bisect();
            self.add_bezier_with_straightness(b0, straightness);
            self.add_bezier_with_straightness(b1, straightness);
        }
    }

    /// Create the flattened segments
    fn create_segments(&mut self) {
        if self.beziers.is_empty() {
            return;
        }
        if !self.flattened.is_empty() {
            return;
        }
        for b in &self.beziers {
            let first = self.flattened.len();
            for pt in b.as_points(self.straightness) {
                if !self.flattened.is_empty() {
                    if vector::distance_sq(&pt, self.flattened.last().unwrap()) > F::epsilon() {
                        self.flattened.push(pt);
                    }
                }
            }
            self.bezier_indices.push((first, self.flattened.len()));
        }
    }
    /// Calculate the distance (squared) from the point to the closest Bezier, and determine
    /// if the point is inside or outside the Beziers if they are a closed path
    pub fn distance_sq_to(&mut self, pt: &[F; 2]) -> (F, bool) {
        self.create_segments();
        let mut min_d2_s = (F::max_value(), 0);
        for (i, (p0, p1)) in self
            .flattened
            .iter()
            .zip(self.flattened.iter().skip(1))
            .enumerate()
        {
            let (distance_sq, inside) = Bezier::pt_distance_sq_from(pt, p0, p1);
            if distance_sq < min_d2_s.0 {
                min_d2_s = (distance_sq, i);
            }
        }
        (min_d2_s.0, false)
    }
}
