use super::{FPType, Fixed, HowIsFixedPoint, UsefulConsts, UsefulInt};

use num_traits::{FromPrimitive, NumCast, ToPrimitive};

impl<T: UsefulInt, const N: usize> FromPrimitive for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn from_u64(n: u64) -> Option<Self> {
        T::from_u64(n).and_then(Self::try_from_t_int)
    }

    fn from_i64(n: i64) -> Option<Self> {
        T::from_i64(n).and_then(Self::try_from_t_int)
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
        let s = self.value >> nb_frac;
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
