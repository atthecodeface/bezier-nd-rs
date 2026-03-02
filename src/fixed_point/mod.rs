//! Fixed point values
mod arith_code;
pub use arith_code::ArithCode;

mod extra_num_traits;
pub use extra_num_traits::{BorrowingSub, CarryingAdd, OverflowingMul};

mod consts;
mod useful;
pub use consts::UsefulConsts;
pub(crate) use useful::{Int, UsefulInt, UsefulUInt};

pub(crate) mod functions;
mod trig;

mod fp_type;
pub use fp_type::{FPType, HowIsFixedPoint};

mod fixed;
pub use fixed::Fixed;

mod binary;
mod unary;

mod nt_casts;
mod nt_consts;
mod nt_float;
mod nt_float_const;
mod nt_num_ops;
mod nt_saturating;
mod nt_signed;

pub trait IsFixed<T, const N: usize> {}
impl<T: UsefulInt, const N: usize> IsFixed<T, N> for Fixed<T, N> where
    FPType<T, N>: HowIsFixedPoint<T>
{
}

mod fp_impls;

#[cfg(test)]
mod test;
