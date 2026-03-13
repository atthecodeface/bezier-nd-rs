use super::{FPType, Fixed, HowIsFixedPoint, UsefulConsts, UsefulInt};

use num_traits::{FromPrimitive, NumCast, ToPrimitive};

impl<T: UsefulInt, const N: usize> FromPrimitive for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn from_u64(n: u64) -> Option<Self> {
        // First generate a T of the value - if this cannot be done, then it won't fit in the integer part of the T
        let Some(v) = T::from_u64(n) else {
            return None;
        };

        // v must be shifted left by NB_FRAC; it must be +ve, since it came from
        // an unsigned value; so if any of the *other* bits are set then the
        // value will not fit
        let nb_frac = <FPType<T, N> as HowIsFixedPoint<T>>::NB_FRAC;
        let nb_int = <FPType<T, N> as HowIsFixedPoint<T>>::NB_INT;
        let nb_max_int_mask = !((!T::ZERO) << nb_int);
        let v_upper_bits = v & !nb_max_int_mask;
        if v_upper_bits.is_zero() {
            Some(Self {
                value: v << nb_frac,
            })
        } else {
            None
        }
    }
    fn from_i64(n: i64) -> Option<Self> {
        // First generate a T of the value - if this cannot be done, then it won't fit in the integer part of the T
        let Some(v) = T::from_i64(n) else {
            return None;
        };

        // v must be shifted left by NB_FRAC; if it is +ve and any of the *other* bits are set then the
        // value will not fit. If it is -ve and any of the *other* bits are clear then the value will not fit
        let nb_frac = <FPType<T, N> as HowIsFixedPoint<T>>::NB_FRAC;
        let nb_int = <FPType<T, N> as HowIsFixedPoint<T>>::NB_INT;
        let nb_max_int_mask = !((!T::ZERO) << nb_int);
        let v_upper_bits = v & !nb_max_int_mask;
        if v_upper_bits.is_zero() || v_upper_bits == !nb_max_int_mask {
            Some(Self {
                value: v << nb_frac,
            })
        } else {
            None
        }
    }
    fn from_f64(n: f64) -> Option<Self> {
        Self::try_from(n).ok()
    }
    fn from_f32(n: f32) -> Option<Self> {
        Self::try_from(n).ok()
    }
}

impl<T: UsefulInt + UsefulConsts, const N: usize> ToPrimitive for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn to_i64(&self) -> Option<i64> {
        let nb_frac = <FPType<T, N> as HowIsFixedPoint<T>>::NB_FRAC;
        let s = self.value >> nb_frac;
        s.to_i64()
    }
    fn to_u64(&self) -> Option<u64> {
        let nb_frac = <FPType<T, N> as HowIsFixedPoint<T>>::NB_FRAC;
        dbg!(self);
        let s = self.value >> nb_frac;
        dbg!(&s);
        s.to_u64()
    }
    fn to_f32(&self) -> Option<f32> {
        Some((*self).into())
    }
    fn to_f64(&self) -> Option<f64> {
        Some((*self).into())
    }
}

impl<T: UsefulInt + UsefulConsts, const N: usize> NumCast for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn from<V: ToPrimitive>(n: V) -> Option<Self> {
        n.to_u64().map(Self::from_u64).flatten()
    }
}
