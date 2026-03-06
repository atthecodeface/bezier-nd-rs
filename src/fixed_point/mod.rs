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

use crate::bignum::{IntN, UIntN};
use num_traits::ConstZero;

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

    // fn unsigned_of_dbl(dbl:Self::Dbl) -> (Self::Unsigned, Self::Unsigned) {
    //     ((dbl >> $nb) as $uns, dbl as $uns)
    // }
}

impl UsefulUInt for UIntN<4> {
    type Dbl = UIntN<8>;
    const NB: usize = 256;
    fn as_dbl(self) -> Self::Dbl {
        (&self).try_into().unwrap()
    }
}

impl UsefulConsts for IntN<4> {
    const E: Self = Self::new(false, [0x_c90f_daa2_2168_c234, 0xc4c6_628b_80dc_1cd1, 0, 0]);
    const LN_2: Self = Self::ZERO;
    const LN_10: Self = Self::ZERO;
    const LOG2_E: Self = Self::ZERO;
    const LOG2_10: Self = Self::ZERO;
    const LOG10_2: Self = Self::ZERO;
    const LOG10_E: Self = Self::ZERO;
    const SQRT_2: Self = Self::ZERO;
    const PI: Self = Self::new(false, [0x_c90f_daa2_2168_c234, 0xc4c6_628b_80dc_1cd1, 0, 0]);
    const FRAC_2_PI: Self = Self::ZERO;
    const FRAC_PI_3: Self = Self::ZERO;
    const FRAC_1_SQRT_PI: Self = Self::ZERO;
}

impl HowIsFixedPoint<IntN<4>> for FPType<IntN<4>, 240> {
    const NB: usize = 256;
    const NB_FRAC: usize = 240;
    const ONE: IntN<4> = IntN::<4>::with_bit_set(240);
}
