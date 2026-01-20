use geo_nd::{vector, Num};

use crate::{BezierEval, BezierSplit, BoxedBezier};

impl<F: 'static + Num + From<f32>, const D: usize> BezierEval<F, [F; D]> for [[F; D]; 2] {
    fn point_at(&self, t: F) -> [F; D] {
        vector::add(vector::scale(self[0], F::ONE - t), &self[1], t)
    }
    fn derivative_at(&self, _t: F) -> [F; D] {
        vector::add(self[1], &self[0], -F::ONE)
    }
    fn endpoints(&self) -> (&[F; D], &[F; D]) {
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
    fn control_point(&self, n: usize) -> &[F; D] {
        &self[n]
    }
}

impl<F: 'static + Num + From<f32>, const D: usize> BezierSplit for [[F; D]; 2] {
    fn split(&self) -> (Self, Self) {
        let m = vector::scale(
            vector::add(self[0], &self[1], 1.0_f32.into()),
            0.5_f32.into(),
        );
        ([self[0], m], [m, self[1]])
    }
}

impl<F: 'static + Num + From<f32>, const D: usize> BezierEval<F, [F; D]> for [[F; D]; 3] {
    fn point_at(&self, t: F) -> [F; D] {
        let two: F = (2.0_f32).into();
        let u = F::ONE - t;
        vector::add(
            vector::add(vector::scale(self[0], u * u), &self[2], t * t),
            &self[1],
            u * t * two,
        )
    }
    fn derivative_at(&self, t: F) -> [F; D] {
        let two: F = (2.0_f32).into();
        let u = F::ONE - t;
        vector::add(
            vector::add(vector::scale(self[0], u * (-two)), &self[2], t * two),
            &self[1],
            (u - t) * two,
        )
    }
    fn endpoints(&self) -> (&[F; D], &[F; D]) {
        (&self[0], &self[1])
    }
    fn is_straight(&self, straightness: F) -> bool {
        let dv = vector::sub(self[2], &self[0], (0.5_f32).into());
        let dv = vector::sub(dv, &self[1], (0.5_f32).into());
        let dc2 = vector::length_sq(&dv);
        dc2 < straightness * straightness
    }
    fn closeness_to_quad(&self) -> F {
        F::ZERO
    }
    fn closeness_to_cubic(&self) -> F {
        F::ZERO
    }
    fn num_control_points(&self) -> usize {
        3
    }
    fn control_point(&self, n: usize) -> &[F; D] {
        &self[n]
    }
}

impl<F: 'static + Num + From<f32>, const D: usize> BoxedBezier<F, [F; D]> for [[F; D]; 2] {
    fn boxed_reduce(&self) -> Option<Box<dyn BoxedBezier<F, [F; D]>>> {
        None
    }
    fn boxed_split(
        &self,
    ) -> Option<(
        Box<dyn BoxedBezier<F, [F; D]>>,
        Box<dyn BoxedBezier<F, [F; D]>>,
    )> {
        let (b0, b1) = <Self as BezierSplit>::split(self);
        Some((Box::new(b0), Box::new(b1)))
    }
}

impl<F: 'static + Num + From<f32>, const D: usize> BezierSplit for [[F; D]; 3] {
    fn split(&self) -> (Self, Self) {
        let c0 = vector::scale(vector::add(self[0], &self[1], F::ONE), (0.5_f32).into());
        let c1 = vector::scale(vector::add(self[1], &self[2], F::ONE), (0.5_f32).into());
        let pm = vector::scale(vector::add(c0, &c1, F::ONE), (0.5_f32).into());
        ([self[0], c0, pm], [pm, c1, self[1]])
    }
}

impl<F: 'static + Num + From<f32>, const D: usize> BoxedBezier<F, [F; D]> for [[F; D]; 3] {
    fn boxed_reduce(&self) -> Option<Box<dyn BoxedBezier<F, [F; D]>>> {
        Some(Box::new([self[0], self[1]]))
    }
    fn boxed_split(
        &self,
    ) -> Option<(
        Box<dyn BoxedBezier<F, [F; D]>>,
        Box<dyn BoxedBezier<F, [F; D]>>,
    )> {
        let (b0, b1) = <Self as BezierSplit>::split(self);
        Some((Box::new(b0), Box::new(b1)))
    }
}
