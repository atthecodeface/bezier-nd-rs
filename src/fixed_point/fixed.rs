use std::sync::atomic::AtomicBool;

use super::UsefulInt;

use super::{ArithCode, FPType, HowIsFixedPoint};

use num_traits::{ConstOne, ConstZero, Zero};

/// A Fixed-point value backed by a signed integer type `T` with `N` bits of
/// fractional part; the value of `N` must be more than zero, and the number of
/// integer bits must be at least one.
#[repr(transparent)]
#[derive(Debug, Default, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Fixed<T: UsefulInt, const N: usize>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    pub(crate) value: T,
}

impl<T: UsefulInt, const N: usize> std::fmt::Display for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        <Self as std::fmt::Debug>::fmt(self, f)
    }
}

impl<T: UsefulInt, const N: usize> std::borrow::Borrow<T> for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn borrow(&self) -> &T {
        &self.value
    }
}

impl<T: UsefulInt, const N: usize> std::borrow::BorrowMut<T> for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn borrow_mut(&mut self) -> &mut T {
        &mut self.value
    }
}

impl<T: UsefulInt, const N: usize> Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    /// Return the underlying value
    ///
    /// Do we want to support Borrow instead? As well?
    pub fn raw(&self) -> T {
        self.value
    }

    /// Take a double and shift it down, using the dropped bits to help round appropriately; return True if it does not overflow,
    /// but always set the value
    ///
    /// Round to nearest, with ties broken *away* from zero (same as f32 or f64.round)
    ///
    /// Away from zero means: for +ve, add 1 if fraction >=1/2, for -ve subtract 1 if fraction >=1/2
    ///
    /// Consider an i16_12, that is 0xd800, so a shift of 8 produces an i8_4, then the i8_4 value should be 0xd8, or -2.5
    ///
    /// 0x2800 represents 2.5
    /// 0x287f represents 2.531005859375 - must round to 2.5 (i.e. no add after shift)
    /// 0x2880 represents 2.53125- must round to 2.5625 (i.e. add 1 after shift)
    /// 0x2880 represents 2.531494140625
    /// 0x2900 represents 2.5625
    /// 0xd800 represents -2.5
    /// 0xd880 represents -2.46875 - must round to -2.5 (i.e. no add after shift)
    /// 0xd881 represents -2.468505859375 - must round to -2.4375 (i.e. add 1 after shift)
    /// 0xd900 represents -2.4375
    ///
    /// Hence if the dropped bits are 0x8000... then add 1 after shift
    fn reduce_double(&mut self, dbl: &T::Dbl, shr: usize) -> bool {
        // Shift left by double size minus the shift right gives an (unsigned) bits we will drop
        // If this is (signed) negative then we are dropping >=1/2; if -itself==itself then it is zero or half
        let bits_to_drop = *dbl & !(!T::Dbl::ZERO << shr);
        let half = T::Dbl::ONE << (shr - 1);
        let bits_to_drop_is_half = bits_to_drop == half;
        let bits_to_drop_is_under_half = (bits_to_drop & half).is_zero();
        let mut add_one_after_shift = !bits_to_drop_is_under_half;
        if *dbl < T::Dbl::ZERO && bits_to_drop_is_half {
            add_one_after_shift = false;
        }

        let mut reduced = *dbl >> shr;
        if add_one_after_shift {
            reduced += T::Dbl::ONE;
        }
        let (value, overflow) = T::of_dbl(reduced);
        self.value = value;
        !overflow
    }

    #[inline(always)]
    pub(crate) fn do_bit_and(&mut self, other: &Self) -> ArithCode {
        self.value &= other.value;
        ArithCode::Ok
    }
    #[inline(always)]
    pub(crate) fn do_bit_or(&mut self, other: &Self) -> ArithCode {
        self.value |= other.value;
        ArithCode::Ok
    }
    #[inline(always)]
    pub(crate) fn do_bit_xor(&mut self, other: &Self) -> ArithCode {
        self.value ^= other.value;
        ArithCode::Ok
    }

    #[inline(always)]
    #[must_use]
    pub(crate) fn do_add(&mut self, other: &Self) -> ArithCode {
        let (r, o) = self.value.overflowing_add(&other.value);
        self.value = r;
        if !o {
            ArithCode::Ok
        } else if self.value < T::ZERO {
            ArithCode::OverflowMin
        } else {
            ArithCode::OverflowMax
        }
    }

    #[inline(always)]
    #[must_use]
    pub(crate) fn do_sub(&mut self, other: &Self) -> ArithCode {
        let (r, o) = self.value.overflowing_sub(&other.value);
        self.value = r;
        if !o {
            ArithCode::Ok
        } else if self.value < T::ZERO {
            ArithCode::OverflowMin
        } else {
            ArithCode::OverflowMax
        }
    }

    // Should probably change to checked mul returning a bool for overflow (and possibly underflow, where result is 0 but neither input is)
    #[inline(always)]
    #[must_use]
    pub(crate) fn do_mul(&mut self, other: &Self) -> ArithCode {
        // If T = x*2^NB_FRAC, then dbl = x*y*2^(NB_FRAC*2), and the result should be x*y*2^NB_FRAC, so shift right by NB_FRAC
        let shr = <FPType<T, N> as HowIsFixedPoint<T>>::NB_FRAC;
        let dbl = self.value.dbl_mult(&other.value);
        if !self.reduce_double(&dbl, shr) {
            if (self.value < T::ZERO) == (other.value < T::ZERO) {
                ArithCode::OverflowMax
            } else {
                ArithCode::OverflowMin
            }
        } else {
            ArithCode::Ok
        }
    }

    // Should probably change to checked div returning a bool for underflow, overflow, div by zero
    #[inline(always)]
    #[must_use]
    pub(crate) fn do_div(&mut self, other: &Self) -> ArithCode {
        if other.value.is_zero() {
            return ArithCode::DivideByZero;
        }
        // If T = x*2^NB_FRAC, then s=x*2^(NB_FRAC+NB), o=y*2^NB_FRAC, r=x/y*2^NB, and the result should be x/y*2^NB_FRAC, so shift right by NB-NB_FRAC = NB_INT
        let shr = <FPType<T, N> as HowIsFixedPoint<T>>::NB_SIGN_AND_INT;

        let s = self.value.as_dbl_upper();
        let o = other.value.as_dbl();
        // r = self/other * 2^(T::NB);
        let r = s / o;

        if !self.reduce_double(&r, shr) {
            if (self.value < T::ZERO) == (other.value < T::ZERO) {
                ArithCode::OverflowMax
            } else {
                ArithCode::OverflowMin
            }
        } else {
            ArithCode::Ok
        }
    }

    // x remainder y = x - (x / y).trunc() * y
    //
    // If x=64 (2^14 in i16_8) and y = 1/256 (1 in i16_8) then x/y = 2^14
    //
    // If x,y are each X or Y*2^N, then x/y = is an integer (i.e. already truncated or rounded in some fashion)
    #[inline(always)]
    pub(crate) fn do_rem(&mut self, other: &Self) -> ArithCode {
        if other.value.is_zero() {
            ArithCode::DivideByZero
        } else {
            // Not conviced this works for negative values
            //
            // And possibly it can just be self.value % other.value...
            let x_div_y_trunc = self.value / other.value;
            self.value -= x_div_y_trunc * other.value;
            ArithCode::Ok
        }
    }
}

#[test]
fn test_thing() {
    let _x = Fixed::<i8, 4>::ONE;
}
