use crate::Num;

/// BezierBuildConstraint is a constraint that can be specified when
/// building a Bezier, to constraint the point at a particular value
/// of parameter 't' to a location, or to constraint the nth derivative
/// of the Bezier at a parameter 't' to a particular value
#[derive(Debug, Clone)]
pub enum BezierBuildConstraint<F: Num, const D: usize> {
    /// Constraint that the Bezier to pass through a particular location at paratmer 't'
    PositionAtT(F, [F; D]),
    /// Constraint that the Bezier to have an nth derivative value at a particular location at paratmer 't'
    DerivativeAtT(F, usize, [F; D]),
}

impl<F: Num, const D: usize> BezierBuildConstraint<F, D> {
    /// Get the depth of derivative of the constraint
    ///
    /// For a constraint of a pure position, this will return 0
    pub fn derivative_depth(&self) -> usize {
        match self {
            BezierBuildConstraint::PositionAtT(_, _) => 0,
            BezierBuildConstraint::DerivativeAtT(_, n, _) => *n,
        }
    }

    /// Get the parameter associated with the constraint
    pub fn at(&self) -> F {
        match self {
            BezierBuildConstraint::PositionAtT(t, _) => *t,
            BezierBuildConstraint::DerivativeAtT(t, _, _) => *t,
        }
    }

    /// Get the position / derivative value associated with the constraint
    pub fn posn(&self) -> &[F; D] {
        match self {
            BezierBuildConstraint::PositionAtT(_, pt) => pt,
            BezierBuildConstraint::DerivativeAtT(_, _, pt) => pt,
        }
    }
}

/// A BezierBuilder is used to create a Bezier with particular constraints,
/// such as positions at certain values of parameter 't', and nth derivatives
/// at other values of 't'
///
/// Once all the constraints are added the Bezier can be built; the degree of
/// the Bezier depends on the number of constraints. Building a Bezier of degree 'N'
/// requires precisely N+1 constraints, and there must not be a constraint for the nth
/// derivative with n greater than the degree.
#[derive(Default, Debug, Clone)]
pub struct BezierBuilder<F: Num, const D: usize> {
    /// The constraints that the Bezier will have to meet after building
    constraints: Vec<BezierBuildConstraint<F, D>>,
}

impl<F: Num, const D: usize> BezierBuilder<F, D> {
    /// Add a constraint that the Bezier passes throught a point at parameter value 't'
    pub fn add_point_at(&mut self, t: F, pt: [F; D]) {
        self.constraints
            .push(BezierBuildConstraint::PositionAtT(t, pt));
    }

    /// Add a constraint that the nth derivative of the Bezier at parameter value 't' has
    /// a particular value
    pub fn add_derivative_at(&mut self, t: F, n: usize, pt: [F; D]) {
        self.constraints
            .push(BezierBuildConstraint::DerivativeAtT(t, n, pt));
    }

    /// Get the number of constraints
    pub fn number_constraints(&self) -> usize {
        self.constraints.len()
    }

    /// Find the largest value of n, for the constraints on the nth derivatives
    ///
    /// If no derivative constraints are present this will be zero
    pub fn deepest_derivative(&self) -> usize {
        let mut deepest = 0;
        for c in &self.constraints {
            deepest = deepest.max(c.derivative_depth());
        }
        deepest
    }

    /// Find the minimum Bezier degree that this requires
    ///
    /// A Bezier of higher degree can be created by elevating a Bezier of the *required*
    /// degree
    pub fn bezier_min_degree(&self) -> usize {
        (self.number_constraints() - 1).max(self.deepest_derivative())
    }

    /// Return an iterator of the build constraints
    pub fn iter(&self) -> impl Iterator<Item = &'_ BezierBuildConstraint<F, D>> {
        self.constraints.iter()
    }
}
