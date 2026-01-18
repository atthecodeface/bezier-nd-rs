//! A simple big-number library providing signed and unsigned integers, and
//! rationals of *fixed* number of 64-bit values
//!
//! It is not highly optimized for performance, and is aimed at being used
//! in algorithms that require numbers that support num_traits::Num

mod int_n;
mod rational_n;
mod uint_n;

pub use int_n::IntN;
pub use rational_n::RationalN;
pub use uint_n::UIntN;
