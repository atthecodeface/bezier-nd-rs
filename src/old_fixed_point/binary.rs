use super::{FPType, Fixed, HowIsFixedPoint, UsefulInt};

use std::ops::*;

macro_rules! binary_op {
    {$trait:ident, $fn:ident, $f:ident, $reason:expr} => {
        impl<T: UsefulInt, const N: usize> $trait<Fixed<T,N>> for Fixed<T, N>
        where
            FPType<T, N>: HowIsFixedPoint<T>,
        {
            type Output = Fixed<T, N>;
            #[track_caller]
            fn $fn (mut self, other:Fixed<T, N>) -> Fixed<T, N> {
                self.$f(&other).check($reason); self
            }
        }
        impl<T: UsefulInt, const N: usize> $trait<&Fixed<T,N>> for Fixed<T, N>
        where
            FPType<T, N>: HowIsFixedPoint<T>,
        {
            type Output = Fixed<T, N>;
            #[track_caller]
            fn $fn (mut self, other:&Fixed<T, N>) -> Fixed<T, N> {
                self.$f(other).check($reason); self
            }
        }
        impl<T: UsefulInt, const N: usize> $trait<Fixed<T,N>> for &Fixed<T, N>
        where
            FPType<T, N>: HowIsFixedPoint<T>,
        {
            type Output = Fixed<T, N>;
            #[track_caller]
            fn $fn (self, other:Fixed<T, N>) -> Fixed<T, N> {
                let mut s = *self; s.$f(&other).check($reason); s
            }
        }
        impl<T: UsefulInt, const N: usize> $trait<&Fixed<T,N>> for &Fixed<T, N>
        where
            FPType<T, N>: HowIsFixedPoint<T>,
        {
            type Output = Fixed<T, N>;
            #[track_caller]
            fn $fn (self, other:&Fixed<T, N>) -> Fixed<T, N> {
                let mut s = *self; s.$f(&other).check($reason); s
            }
        }
    }
}

macro_rules! binary_assign_op {
    {$trait:ident, $fn:ident, $f:ident, $reason:expr} => {
    impl<T: UsefulInt, const N: usize> $trait<Fixed<T,N>> for Fixed<T, N>
    where
        FPType<T, N>: HowIsFixedPoint<T>,
    {
            #[track_caller]
            fn $fn (&mut self, other:Fixed<T, N>) {
                self.$f(&other).check($reason);
            }
        }
        impl<T: UsefulInt, const N: usize> $trait<&Fixed<T,N>> for Fixed<T, N>
        where
            FPType<T, N>: HowIsFixedPoint<T>,
        {
            #[track_caller]
            fn $fn (&mut self, other:&Fixed<T, N>) {
                self.$f(other).check($reason);
            }
        }
    }
}

binary_op! {Add, add, do_add, "addition"}
binary_assign_op! {AddAssign, add_assign, do_add, "addition"}
binary_op! {Sub, sub, do_sub, "subtraction"}
binary_assign_op! {SubAssign, sub_assign, do_sub, "subtraction"}
binary_op! {Mul, mul, do_mul, "multiplcation"}
binary_assign_op! {MulAssign, mul_assign, do_mul, "multiplcation"}
binary_op! {Div, div, do_div, "division"}
binary_assign_op! {DivAssign, div_assign, do_div, "division"}
binary_op! {Rem, rem, do_rem, "remainder"}
binary_assign_op! {RemAssign, rem_assign, do_rem, "remainder"}
binary_op! {BitAnd, bitand, do_bit_and, ""}
binary_assign_op! {BitAndAssign, bitand_assign, do_bit_and, ""}
binary_op! {BitOr, bitor, do_bit_or, ""}
binary_assign_op! {BitOrAssign, bitor_assign, do_bit_or, ""}
binary_op! {BitXor, bitxor, do_bit_xor, ""}
binary_assign_op! {BitXorAssign, bitxor_assign, do_bit_xor, ""}
