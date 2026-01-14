//a Imports
use geo_nd::Float;

//a Bezier
pub enum BezierBuilder<F: Float, const D: usize> {
    PositionAtT(F, [F; D]),
    DerivativeAtT(F, usize, [F; D]),
}
