//a Documentation
#![warn(missing_docs)]
#![doc = include_str!("lib_doc.md")]

/*a Imports
*/
mod error;
pub use error::BezierError;
mod builder;
/// Constants used for Bezier curves and their elevation/reduction
pub mod constants;
pub mod lazy_constants;
pub(crate) mod utils;

mod traits;
pub use traits::{Num, NumOps};

mod approximation;
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
pub use bezier_iter::{BezierLineTIter, BezierPointTIter, BezierSplitTIter};
pub use traits::{
    BasicBezier, BezierConstruct, BezierElevate, BezierEval, BezierFlatIterator,
    BezierIterationType, BezierMetric, BezierOps, BezierReduce, BezierReduction, BezierSplit,
    BoxedBezier,
};

pub use builder::BezierBuilder;
pub use constants::*;
pub use polynomial::{PolyFindRoots, PolyNewtonRaphson, Polynomial};

// mod distance;
// pub use distance::BezierDistance2D;

pub use implementations::bezier_nd::BezierND;
