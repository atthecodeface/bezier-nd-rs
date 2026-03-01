use super::{FPType, Fixed, HowIsFixedPoint, UsefulConsts, UsefulInt};

use num_traits::FloatConst;

use std::ops::*;

macro_rules! binary_op {
    {$trait:ident, $fn:ident, $f:ident} => {
        impl<T: UsefulInt, const N: usize> $trait<Fixed<T,N>> for Fixed<T, N>
        where
            FPType<T, N>: HowIsFixedPoint<T>,
        {
            type Output = Fixed<T, N>;
            fn $fn (mut self, other:Fixed<T, N>) -> Fixed<T, N> {
                self.$f(&other); self
            }
        }
        impl<T: UsefulInt, const N: usize> $trait<&Fixed<T,N>> for Fixed<T, N>
        where
            FPType<T, N>: HowIsFixedPoint<T>,
        {
            type Output = Fixed<T, N>;
            fn $fn (mut self, other:&Fixed<T, N>) -> Fixed<T, N> {
                self.$f(other); self
            }
        }
        impl<T: UsefulInt, const N: usize> $trait<Fixed<T,N>> for &Fixed<T, N>
        where
            FPType<T, N>: HowIsFixedPoint<T>,
        {
            type Output = Fixed<T, N>;
            fn $fn (self, other:Fixed<T, N>) -> Fixed<T, N> {
                let mut s = *self; s.$f(&other); s
            }
        }
        impl<T: UsefulInt, const N: usize> $trait<&Fixed<T,N>> for &Fixed<T, N>
        where
            FPType<T, N>: HowIsFixedPoint<T>,
        {
            type Output = Fixed<T, N>;
            fn $fn (self, other:&Fixed<T, N>) -> Fixed<T, N> {
                let mut s = *self; s.$f(&other); s
            }
        }
    }
}

macro_rules! binary_assign_op {
    {$trait:ident, $fn:ident, $f:ident} => {
    impl<T: UsefulInt, const N: usize> $trait<Fixed<T,N>> for Fixed<T, N>
    where
        FPType<T, N>: HowIsFixedPoint<T>,
    {
            fn $fn (&mut self, other:Fixed<T, N>) {
                self.$f(&other);
            }
        }
        impl<T: UsefulInt, const N: usize> $trait<&Fixed<T,N>> for Fixed<T, N>
        where
            FPType<T, N>: HowIsFixedPoint<T>,
        {
            fn $fn (&mut self, other:&Fixed<T, N>) {
                self.$f(other);
            }
        }
    }
}

binary_op! {Add, add, do_add}
binary_assign_op! {AddAssign, add_assign, do_add}
binary_op! {Sub, sub, do_sub}
binary_assign_op! {SubAssign, sub_assign, do_sub}
binary_op! {Mul, mul, do_mul}
binary_assign_op! {MulAssign, mul_assign, do_mul}
binary_op! {Div, div, do_div}
binary_assign_op! {DivAssign, div_assign, do_div}
binary_op! {Rem, rem, do_rem}
binary_assign_op! {RemAssign, rem_assign, do_rem}
binary_op! {BitAnd, bitand, do_bit_and}
binary_assign_op! {BitAndAssign, bitand_assign, do_bit_and}
binary_op! {BitOr, bitor, do_bit_or}
binary_assign_op! {BitOrAssign, bitor_assign, do_bit_or}
binary_op! {BitXor, bitxor, do_bit_xor}
binary_assign_op! {BitXorAssign, bitxor_assign, do_bit_xor}
