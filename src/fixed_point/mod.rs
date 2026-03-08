//! Fixed point values
mod arith_code;
pub use arith_code::ArithCode;

mod consts;
mod useful;
pub use consts::{u64_4, UsefulConsts};
pub use useful::{Int, UsefulInt, UsefulUInt};

pub(crate) mod functions;
mod trig;

mod fp_type;
pub use fp_type::{FPType, HowIsFixedPoint};

mod fixed;
pub use fixed::Fixed;

mod binary;
mod unary;

mod nt_casts;
mod nt_checked;
mod nt_consts;
mod nt_float;
mod nt_float_const;
mod nt_num_ops;
mod nt_overflowing;
mod nt_saturating;
mod nt_signed;
mod nt_wrapping;

mod nt_num;

/// A marker trait that is implemented for `Fixed<T,N>` for any type T providing a fixed point value with N fraction bits
///
/// This is implemented automatically when the trait `HowIsFixedPoint<T>` is implemented for `FPType<T,N>`.
///
/// Quite probably the [UsefulConsts] should be somehow part of this? But they are for T...
///
/// This should probably be sealed as external types *cannot* implement it?
pub trait IsFixed<T, const N: usize> {}
impl<T: UsefulInt, const N: usize> IsFixed<T, N> for Fixed<T, N> where
    FPType<T, N>: HowIsFixedPoint<T>
{
}

mod fp_impls;

#[cfg(test)]
mod test;

use crate::bignum::{IntN, UIntN};
use num_traits::ConstZero;

// It would be nice to punt this into bignum; then it proves that Fixed can run
// on user types. It is probably better that way round than Fixed running only
// on types that it knows about
impl UsefulInt for IntN<4> {
    type Unsigned = UIntN<4>;
    type Dbl = IntN<8>;
    type DblUnsigned = UIntN<8>;
    const NB: usize = 256;
    fn as_dbl(self) -> Self::Dbl {
        let v: UIntN<8> = (&self.value).try_into().unwrap();
        v.into()
    }
    fn as_dbl_upper(self) -> Self::Dbl {
        self.as_dbl() << Self::NB
    }
    fn as_dbl_unsigned(self) -> Self::DblUnsigned {
        (&self.value).try_into().unwrap()
    }
    fn unsigned(self) -> Self::Unsigned {
        self.value
    }
    fn of_unsigned(v: Self::Unsigned) -> Self {
        v.into()
    }
    fn dbl_mult(self, other: &Self) -> Self::Dbl {
        self.as_dbl() * other.as_dbl()
    }
    fn of_dbl(dbl: Self::Dbl) -> (Self, bool) {
        let value = dbl.value & UIntN::<8>::mask(Self::NB as u32);
        let overflow = value != dbl.value;
        let mut s = Self::default();
        s.value = (&value).try_into().unwrap();
        (s, overflow)
    }
    fn min_value() -> Self {
        Self {
            is_neg: true,
            value: !UIntN::<4>::ZERO,
        }
    }
    fn max_value() -> Self {
        Self {
            is_neg: false,
            value: !UIntN::<4>::ZERO,
        }
    }
}

impl UsefulUInt for UIntN<4> {
    type Dbl = UIntN<8>;
    const NB: usize = 256;
    fn as_dbl(self) -> Self::Dbl {
        (&self).try_into().unwrap()
    }
}

impl UsefulConsts for IntN<4> {
    // These constants must be aligned to the top bit; the number of fractional bits is taken into account by the 'Fixed' type
    const E: Self = Self::new(false, consts::u64_4::E);
    const LN_2: Self = Self::new(false, consts::u64_4::LN_2);
    const LN_10: Self = Self::new(false, consts::u64_4::LN_10);
    const LOG2_E: Self = Self::new(false, consts::u64_4::LOG2_E);
    const LOG2_10: Self = Self::new(false, consts::u64_4::LOG2_10);
    const LOG10_E: Self = Self::new(false, consts::u64_4::LOG10_E);
    const LOG10_2: Self = Self::new(false, consts::u64_4::LOG10_2);
    const SQRT_2: Self = Self::new(false, consts::u64_4::SQRT_2);
    const PI: Self = Self::new(false, consts::u64_4::PI);
    const FRAC_2_PI: Self = Self::new(false, consts::u64_4::FRAC_2_PI);
    const FRAC_PI_3: Self = Self::new(false, consts::u64_4::FRAC_PI_3);
    const FRAC_1_SQRT_PI: Self = Self::new(false, consts::u64_4::FRAC_1_SQRT_PI);
}

impl HowIsFixedPoint<IntN<4>> for FPType<IntN<4>, 240> {
    const NB: usize = 256;
    const NB_FRAC: usize = 240;
    const NB_INT: usize = Self::NB - Self::NB_FRAC;
    const ONE: IntN<4> = IntN::<4>::with_bit_set(240);
    const DEDICATED_SIGN: bool = true;
}
