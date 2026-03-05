use super::IsFixed;

trait TestKind<I, const N: usize>:
    num_traits::Num
    + num_traits::ConstOne
    + num_traits::ConstZero
    + std::fmt::Debug
    + num_traits::FromPrimitive
    + PartialOrd
    + Copy
    + std::borrow::Borrow<I>
    + std::borrow::BorrowMut<I>
    + IsFixed<I, N>
{
}
impl<T, I, const N: usize> TestKind<I, N> for T where
    T: num_traits::Num
        + num_traits::ConstOne
        + num_traits::ConstZero
        + std::fmt::Debug
        + num_traits::FromPrimitive
        + PartialOrd
        + Copy
        + std::borrow::Borrow<I>
        + std::borrow::BorrowMut<I>
        + IsFixed<I, N>
{
}

trait BackingKind:
    num_traits::Num
    + num_traits::ConstOne
    + num_traits::ConstZero
    + std::fmt::Debug
    + num_traits::FromPrimitive
    + num_traits::ToPrimitive
    + PartialOrd
    + Copy
{
}

impl<T> BackingKind for T where
    T: num_traits::Num
        + num_traits::ConstOne
        + num_traits::ConstZero
        + std::fmt::Debug
        + num_traits::FromPrimitive
        + num_traits::ToPrimitive
        + PartialOrd
        + Copy
{
}

pub(self) mod calculate_constants;
pub(self) mod constants;
mod test_i8;
pub(self) mod utils;
