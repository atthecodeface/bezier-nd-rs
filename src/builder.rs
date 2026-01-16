use std::{fmt::Debug, ops::Deref};

//a Imports
use geo_nd::Float;

//a Bezier
#[derive(Debug, Clone)]
pub enum BezierBuildConstraint<F: Float, const D: usize> {
    PositionAtT(F, [F; D]),
    DerivativeAtT(F, usize, [F; D]),
}

impl<F: Float, const D: usize> BezierBuildConstraint<F, D> {
    pub fn derivative_depth(&self) -> usize {
        match self {
            BezierBuildConstraint::PositionAtT(_, _) => 0,
            BezierBuildConstraint::DerivativeAtT(_, n, _) => *n,
        }
    }
    pub fn at(&self) -> F {
        match self {
            BezierBuildConstraint::PositionAtT(t, _) => *t,
            BezierBuildConstraint::DerivativeAtT(t, _, _) => *t,
        }
    }
    pub fn posn(&self) -> &[F; D] {
        match self {
            BezierBuildConstraint::PositionAtT(_, pt) => pt,
            BezierBuildConstraint::DerivativeAtT(_, _, pt) => pt,
        }
    }
}

#[derive(Default, Debug, Clone)]
pub struct BezierBuilder<F: Float, const D: usize> {
    constraints: Vec<BezierBuildConstraint<F, D>>,
}

impl<F: Float, const D: usize> BezierBuilder<F, D> {
    pub fn add_point_at(&mut self, t: F, pt: [F; D]) {
        self.constraints
            .push(BezierBuildConstraint::PositionAtT(t, pt));
    }
    pub fn number_constraints(&self) -> usize {
        self.constraints.len()
    }
    pub fn deepest_derivative(&self) -> usize {
        let mut deepest = 0;
        for c in &self.constraints {
            deepest = deepest.max(c.derivative_depth());
        }
        deepest
    }
    pub fn bezier_min_degree(&self) -> usize {
        (self.number_constraints() - 1).max(self.deepest_derivative())
    }
    pub fn iter(&self) -> impl Iterator<Item = &'_ BezierBuildConstraint<F, D>> {
        self.constraints.iter()
    }
}
