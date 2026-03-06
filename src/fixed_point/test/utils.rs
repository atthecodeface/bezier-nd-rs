use super::{BackingKind, TestKind};

#[track_caller]
pub fn raw_u64<T, I, const N: usize>(t: &T) -> u64
where
    T: TestKind<I, N>,
    I: BackingKind,
{
    t.borrow().to_u64().unwrap()
}

pub fn set_raw<T, I, const N: usize>(t: &mut T, value: u64)
where
    T: TestKind<I, N>,
    I: BackingKind,
{
    *(t.borrow_mut()) = I::from_u64(value).unwrap();
}

pub fn int_bits<T, I, const N: usize>(_t: &T) -> usize
where
    T: TestKind<I, N>,
    I: BackingKind,
{
    std::mem::size_of::<I>() * 8 - N - 1
}

pub fn nearly_equal<T, I, const N: usize>(t: T, e: T, frac_value: u64) -> bool
where
    T: TestKind<I, N> + num_traits::Float,
    I: BackingKind,
{
    let diff = raw_u64(&((t - e).abs()));
    diff <= frac_value
}
