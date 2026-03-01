use super::{FPType, Fixed, HowIsFixedPoint, UsefulInt};

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
        let nb_max_int_mask = <FPType<T, N> as HowIsFixedPoint<T>>::MAX_INT_MASK;
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
        let nb_max_int_mask = <FPType<T, N> as HowIsFixedPoint<T>>::MAX_INT_MASK;
        let v_upper_bits = v & !nb_max_int_mask;
        if v_upper_bits.is_zero() || v_upper_bits == !nb_max_int_mask {
            Some(Self {
                value: v << nb_frac,
            })
        } else {
            None
        }
    }
}

impl<T: UsefulInt, const N: usize> ToPrimitive for Fixed<T, N>
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
        let s = self.value >> nb_frac;
        s.to_u64()
    }
}

impl<T: UsefulInt, const N: usize> NumCast for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn from<V: ToPrimitive>(n: V) -> Option<Self> {
        n.to_u64().map(Self::from_u64).flatten()
    }
}
