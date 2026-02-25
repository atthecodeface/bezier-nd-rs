use super::{FixedPoint, SignedRaw3232};

#[allow(non_camel_case_types)]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
/// Fixed point i32 with 16 bits of fractional value
pub struct FixedPoint_i32_16 {
    value: i32,
}

impl FixedPoint_i32_16 {
    pub(crate) const fn of_i8(v: i8) -> Self {
        Self {
            value: (v as i32) << Self::FRAC_BITS,
        }
    }
}

impl std::convert::From<(i16, u16)> for FixedPoint_i32_16 {
    fn from((i, f): (i16, u16)) -> Self {
        Self {
            value: ((i as i32) << 16) | (f as i32),
        }
    }
}

impl std::convert::TryFrom<SignedRaw3232> for FixedPoint_i32_16 {
    type Error = ();
    fn try_from(value: SignedRaw3232) -> Result<Self, ()> {
        let i = value.int();
        if i <= (i16::MAX as u32) {
            let i = i as i16;
            let f = (value.frac() >> 16) as u16;
            if value.is_neg() {
                Ok((-i, f).into())
            } else {
                Ok((-i, f).into())
            }
        } else {
            Err(())
        }
    }
}

impl std::convert::From<FixedPoint_i32_16> for SignedRaw3232 {
    fn from(value: FixedPoint_i32_16) -> Self {
        (value.value >> 16, (value.value & 0xffff) as u32).into()
    }
}

// Arithemtic operations used by macros
impl FixedPoint_i32_16 {
    #[inline(always)]
    pub(crate) fn do_add(&mut self, other: &Self) {
        self.value += other.value;
    }
    #[inline(always)]
    pub(crate) fn do_sub(&mut self, other: &Self) {
        self.value -= other.value;
    }
    #[inline(always)]
    pub(crate) fn do_bit_and(&mut self, other: &Self) {
        self.value &= other.value;
    }
    #[inline(always)]
    pub(crate) fn do_bit_or(&mut self, other: &Self) {
        self.value |= other.value;
    }
    #[inline(always)]
    pub(crate) fn do_bit_xor(&mut self, other: &Self) {
        self.value ^= other.value;
    }
    #[inline(always)]
    pub(crate) fn do_mul(&mut self, other: &Self) {
        let s: SignedRaw3232 = (*self).into();
        let o: SignedRaw3232 = (*other).into();
        *self = (s * o).try_into().unwrap();
    }
    #[inline(always)]
    pub(crate) fn do_div(&mut self, other: &Self) {
        let s: SignedRaw3232 = (*self).into();
        let o: SignedRaw3232 = (*other).into();
        *self = (s / o).try_into().unwrap();
    }
    #[inline(always)]
    pub(crate) fn do_rem(&mut self, other: &Self) {
        let s: SignedRaw3232 = (*self).into();
        let o: SignedRaw3232 = (*other).into();
        *self = (s % o).try_into().unwrap();
    }
}

impl std::fmt::Display for FixedPoint_i32_16 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.value.fmt(f)
    }
}

impl FixedPoint for FixedPoint_i32_16 {
    const FRAC_BITS: usize = 16;
    const INT_BITS: usize = 16;
    const SIGNED: bool = true;
    type INT = i32;
    type FRAC = u32;
    type DBLINT = i32;
    type DBLFRAC = u32;
    type RAW = i32;

    #[inline(always)]
    fn as_raw(&self) -> Self::RAW {
        self.value
    }

    #[inline(always)]
    fn of_raw(value: &Self::RAW) -> Self {
        Self {
            value: (*value as i32),
        }
    }

    fn inter_int(&self) -> Self::DBLINT {
        self.value >> Self::FRAC_BITS
    }
    fn inter_frac(&self) -> Self::DBLFRAC {
        (self.value as Self::DBLFRAC) & (Self::FRAC_MASK as Self::DBLFRAC)
    }

    fn of_inter(value: &(Self::DBLINT, Self::DBLFRAC)) -> Option<Self> {
        let f = value.1;
        let i = value.0 + (f >> Self::FRAC_BITS) as Self::DBLINT;
        let result = ((i as Self::RAW) << Self::INT_BITS)
            + ((f as Self::RAW) & (Self::FRAC_MASK as Self::RAW));
        let result = Self { value: result };
        if result.inter_int() != i {
            None
        } else {
            Some(result)
        }
    }
}
