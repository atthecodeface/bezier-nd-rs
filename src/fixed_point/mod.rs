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
    // These constants must be aligned to the top bit; the number of fractional bits is taken into account by the 'Fixed' type
    const E: Self = Self::new(
        false,
        [
            0xadf85458a2bb4a9a,
            0xafdc5620273d3cf1,
            0xd8b9c583ce2d3695,
            0xa9e13641146433f4,
        ],
    );
    const LN_2: Self = Self::new(
        false,
        [
            0xb17217f7d1cf79ab,
            0xc9e3b39803f2f6af,
            0x40f343267298b62d,
            0x8a0d175b8baafa2b,
        ],
    );
    const LN_10: Self = Self::new(
        false,
        [
            0x935d8dddaaa8ac16,
            0xea56d62b82d30a28,
            0xe28fecf9da5df90e,
            0x83c61e8201f02d72,
        ],
    );
    const LOG2_E: Self = Self::new(
        false,
        [
            0xb8aa3b295c17f0bb,
            0xbe87fed0691d3e88,
            0xeb577aa8dd695a58,
            0x8b25166cd1a13247,
        ],
    );
    const LOG2_10: Self = Self::new(
        false,
        [
            0xd49a784bcd1b8afe,
            0x492bf6ff4dafdb4c,
            0xd96c55fe37b3ad4e,
            0x91b6ac8082e7859d,
        ],
    );
    const LOG10_E: Self = Self::new(
        false,
        [
            0xde5bd8a937287195,
            0x355baaafad33dc32,
            0x3ee3460245c9a202,
            0x3a3f2d44f78ea53c,
        ],
    );
    const LOG10_2: Self = Self::new(
        false,
        [
            0x9a209a84fbcff798,
            0x8f8959ac0b7c9178,
            0x26ad30c543d1f349,
            0x8a5e6f26b7cc63cb,
        ],
    );
    const SQRT_2: Self = Self::new(
        false,
        [
            0xb504f333f9de6484,
            0x597d89b3754abe9f,
            0x1d6f60ba893ba84c,
            0xed17ac8583339915,
        ],
    );
    const PI: Self = Self::new(
        false,
        [
            0xc90fdaa22168c234,
            0xc4c6628b80dc1cd1,
            0x29024e088a67cc74,
            0x20bbea63b139b1a,
        ],
    );
    const FRAC_2_PI: Self = Self::new(
        false,
        [
            0xa2f9836e4e441529,
            0xfc2757d1f534ddc0,
            0xdb6295993c439041,
            0xfe5163abdebbc568,
        ],
    );
    const FRAC_PI_3: Self = Self::new(
        false,
        [
            0x860a91c16b9b2c23,
            0x2dd99707ab3d688b,
            0x70ac3405b19a884d,
            0x56b27f197cb7bcbc,
        ],
    );
    const FRAC_1_SQRT_PI: Self = Self::new(
        false,
        [
            0x906eba8214db688d,
            0x71d48a7f6bfec344,
            0x1409a0ebac3e7517,
            0x39a15830cce620b0,
        ],
    );
}

impl HowIsFixedPoint<IntN<4>> for FPType<IntN<4>, 240> {
    const NB: usize = 256;
    const NB_FRAC: usize = 240;
    const NB_INT: usize = Self::NB - Self::NB_FRAC;
    const ONE: IntN<4> = IntN::<4>::with_bit_set(240);
}
