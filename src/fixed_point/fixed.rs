use crate::fixed_point::functions::i32_28::HALF_PI_U32_28;

use super::{UsefulConsts, UsefulInt};

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
    pub fn raw(&self) -> &T {
        &self.value
    }

    /// Return the underlying value
    ///
    /// Do we want to support Borrow instead? As well?
    pub fn raw_mut(&mut self) -> &mut T {
        &mut self.value
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

    /// Calculate which quadrant, and what fraction within the octant to use for sin/cos
    ///
    /// The quadrants are PI/2 around each of the four axes; quadrant 0 is thus
    /// from 'self' being [-PI/4, +PI/4) The fraction returned is [-1, 1). The
    /// fraction returned is always in T as a complete fraction.
    ///
    /// If T is twos complement:
    ///  Since self is angle * 2^N, and PI is stored in [UsefulConsts]
    ///  as PI * 2^(T::NB-3), if we calculate self * 2^(T::NB) / (PI * 2^(T::NB-3))
    ///  we will get angle/PI * 2^(N+3)
    ///
    ///  The octant is therefore in bits (N+1)..+3, and the fraction of PI/2 is in bits 1..(N)
    pub fn quadrant_angle(&self) -> (u8, T)
    where
        T: UsefulConsts,
    {
        use num_traits::ToPrimitive;
        let s = self.value.as_dbl_upper();
        // n_pi is actually PI*2^(T::NB-3) for T being 2s complement, PI*2^(T::NB-2) for T with a dedicated sign bit
        let n_pi = T::PI.as_dbl();
        let r = s / n_pi;
        // Select the correct 3 bits as the octant
        let octant = (r >> (N + 1)) & !(!T::Dbl::ZERO << 3_usize);
        let octant = octant.to_u8().unwrap();
        let quadrant = ((octant + 1) & 6) >> 1;
        // Select the angle - note this will have the top bit (of T) set for odd
        // octants, so it is a negative value, which is what we want
        let angle = (r << (T::NB - (N + 2))) & !(!T::Dbl::ZERO << T::NB);
        (quadrant, T::of_dbl(angle).0)
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
    /// Calculate sin and cos of an angle in the range -PI/2 to +PI/2
    ///
    /// Should clamp in range
    pub fn sincos_first_quad(angle: i32) -> (i32, i32) {
        use super::functions::i32_28::{
            ATAN_ANGLES_I32_28, COS_SCALE_U32_28, HALF_PI_U32_28, NEG_POW2_I32_28,
        };
        let (c, s) = super::functions::apply_rotation_table::<_, 1>(
            angle,
            (COS_SCALE_U32_28 as i32, 0),
            ATAN_ANGLES_I32_28,
            NEG_POW2_I32_28,
        )
        .1;
        (s, c)
    }

    use num_traits::{FloatConst, FromPrimitive};
    for n in 0..100 {
        let angle_f = (n as f32) * std::f32::consts::PI / 50.0;
        let mut angle = Fixed::<i32, 16>::from_f32(angle_f).unwrap();
        *angle.raw_mut() = (angle_f * 65536.0) as i32;
        eprintln!(
            "Angle: {angle_f} {angle} {} sin:{} cos:{}",
            angle.value as f32 / (65536.0),
            angle_f.sin(),
            angle_f.cos()
        );
        let (quadrant, subangle) = angle.quadrant_angle();
        eprintln!(
            "{quadrant} {} {}",
            subangle as f32 / 65536.0,
            subangle as f32 / (65536.0 * 32768.0)
        );
        *angle.raw_mut() = subangle;
        let a = angle.value.as_dbl() * (HALF_PI_U32_28 as i64);
        let sincos = sincos_first_quad((a >> 32) as i32);
        let (sin, cos) = match quadrant {
            0 => (sincos.0, sincos.1),
            1 => (sincos.1, -sincos.0),
            2 => (-sincos.0, -sincos.1),
            _ => (-sincos.1, sincos.0),
        };
        eprintln!(
            "{} {} {}",
            a as f32 / (1_u64 << 44) as f32,
            sin as f32 / (1_u64 << 28) as f32,
            cos as f32 / (1_u64 << 28) as f32,
        );
    }
    assert!(false, "Force failure");
}
