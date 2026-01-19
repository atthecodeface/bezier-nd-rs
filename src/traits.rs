use geo_nd::{vector, Num};

pub trait DynBezier<F: Num> {
    type Point: Clone;
    fn point_at(&self, t: F) -> Self::Point;
    fn derivative_at(&self, t: F) -> Self::Point;
    fn endpoints(&self) -> (&Self::Point, &Self::Point);
    fn is_straight(&self, straightness: F) -> bool;
    fn closeness_to_quad(&self) -> F;
    fn closeness_to_cubic(&self) -> F;
    fn num_control_points(&self) -> usize;
    fn control_point(&self, n: usize) -> &Self::Point;
}
pub trait BezierSplit: Sized {
    fn bisect(&self) -> (Self, Self);
}

impl<F: Num + From<f32>, const D: usize> DynBezier<F> for [[F; D]; 2] {
    type Point = [F; D];
    fn point_at(&self, t: F) -> [F; D] {
        vector::add(vector::scale(self[0], F::ONE - t), &self[1], t)
    }
    fn derivative_at(&self, _t: F) -> [F; D] {
        vector::add(self[1], &self[0], -F::ONE)
    }
    fn endpoints(&self) -> (&Self::Point, &Self::Point) {
        (&self[0], &self[1])
    }
    fn is_straight(&self, _straightness: F) -> bool {
        true
    }
    fn closeness_to_quad(&self) -> F {
        F::ZERO
    }
    fn closeness_to_cubic(&self) -> F {
        F::ZERO
    }
    fn num_control_points(&self) -> usize {
        2
    }
    fn control_point(&self, n: usize) -> &Self::Point {
        &self[n]
    }
}

impl<F: Num + From<f32>, const D: usize> BezierSplit for [[F; D]; 2] {
    fn bisect(&self) -> (Self, Self) {
        let m = vector::scale(
            vector::add(self[0], &self[1], 1.0_f32.into()),
            0.5_f32.into(),
        );
        ([self[0], m], [m, self[1]])
    }
}
