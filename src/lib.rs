#![warn(missing_docs)]
#![doc = include_str!("lib_doc.md")]

/*a Imports
*/
mod error;
pub use error::BezierError;
mod builder;

/// Constants used for Bezier curves and their elevation/reduction
pub mod constants;
mod constants_table;

pub(crate) mod utils;

mod traits;

mod approximation;

pub(crate) mod polynomial;

mod bezier_iter;
mod implementations;

/// Useful functions for generating Bernstein coefficients and operating on Berstein Bezier curves
pub mod bernstein_fns;

pub mod metrics;

/*a Exports
*/
pub use approximation::Approximation;
pub(crate) use bezier_iter::{BezierLineTIter, BezierPointTIter, BezierSplitTIter};
pub use constants_table::{
    ConstantsTable, SharedConstantsTable, StaticConstantsTable, CONSTANTS_TABLE,
};
pub use traits::Num;
pub use traits::{
    BezierConstruct, BezierElevate, BezierEval, BezierFlatIterator, BezierIterationType, BezierMap,
    BezierMetric, BezierOps, BezierReduce, BezierReduction, BezierSplit,
};

pub use builder::BezierBuilder;
pub use constants::*;
pub use polynomial::{PolyFindRoots, PolyNewtonRaphson, Polynomial};

// mod distance;
// pub use distance::BezierDistance2D;

pub use implementations::bezier_nd::BezierND;

impl<T: fixed_pt::UsefulInt + fixed_pt::UsefulConsts, const N: usize> crate::Num
    for fixed_pt::Fixed<T, N>
where
    fixed_pt::IsAFixedType<T, N>: fixed_pt::HowIsFixed<T>,
{
    /// Return true if the divisor is too close to zero
    /// to provide a useful result
    fn is_unreliable_divisor(self) -> bool {
        self.raw().is_zero()
    }

    /// Return an estimate of the square root to within a precision
    fn sqrt_est(self) -> Self {
        // sqrt with a u32_24 input produces a u32_28 output
        // ; a u64_0 (pure integer) produces a u64_32 output
        let nb = <fixed_pt::IsAFixedType<T, N> as fixed_pt::HowIsFixed<T>>::NB_FRAC;
        if (nb & 1) == 0 {
            let x_uns_half_frac = fixed_pt::functions::sqrt(self.raw().unsigned(), false);
            Self::of_raw(T::of_unsigned(x_uns_half_frac) >> (nb / 2))
        } else {
            let x_uns_half_frac = fixed_pt::functions::sqrt(self.raw().unsigned(), true);
            Self::of_raw(T::of_unsigned(x_uns_half_frac) >> (nb / 2))
        }
    }

    /// Return true if the divisor is too close to zero
    /// to provide a useful result
    fn cbrt_est(self) -> Self {
        crate::utils::cbrt_est::<_, 5>(self)
    }
}

#[test]
fn test_me() {
    let _x = BezierND::<fixed_pt::Fixed<i32, 16>, 2, 1>::default();
}
