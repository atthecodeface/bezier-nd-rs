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

mod nt_binary;
mod nt_casts;
mod nt_consts;
mod nt_float;
mod nt_float_const;
mod nt_num_ops;
mod nt_saturating;
mod nt_signed;
mod nt_unary;
