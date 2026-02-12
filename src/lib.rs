//a Documentation
#![warn(missing_docs)]
#![doc = include_str!("lib_doc.md")]

/*a Imports
*/
mod builder;
pub(crate) mod constants;
pub(crate) mod utils;

mod traits;
pub use traits::{Float, Num};

mod approximation;
mod distance;
pub(crate) mod polynomial;

mod bezier_iter;

mod implementations;

/// Useful functions for generating Bernstein coefficients and operating on Berstein Bezier curves
pub mod bernstein_fns;

pub mod bignum;

pub mod metrics;

/*a Exports
*/
pub use approximation::Approximation;
pub use bezier_iter::{
    BezierLineIter, BezierLineTIter, BezierPointIter, BezierPointTIter, BezierQuadIter,
    BezierSplitIter, BezierSplitTIter,
};
pub use traits::{
    BasicBezier, BezierConstruct, BezierElevate, BezierEval, BezierFlatIterator, BezierMetric,
    BezierOps, BezierReduce, BezierReduction, BezierSection, BezierSplit, BoxedBezier,
};

pub use builder::BezierBuilder;
pub use constants::*;
pub use distance::BezierDistance2D;
pub use polynomial::{PolyFindRoots, PolyNewtonRaphson, Polynomial};

pub use implementations::bezier_nd::BezierND;
