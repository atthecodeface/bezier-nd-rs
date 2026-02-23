use super::FixedPoint;
use super::{FixedPoint_i32_16, SignedRaw3232, UnsignedRaw3232};

macro_rules! binary_op {
    {$t:ty, $trait:ident, $fn:ident, $f:ident} => {
        impl $trait < $t > for $t {
            type Output = $t;
            fn $fn (mut self, other:$t) -> $t {
                self.$f(&other); self
            }
        }
        impl $trait < &$t > for $t {
            type Output = $t;
            fn $fn (mut self, other:&$t) -> $t {
                self.$f(other); self
            }
        }
        impl $trait < $t > for &$t {
            type Output = $t;
            fn $fn (self, other:$t) -> $t {
                let mut s = *self; s.$f(&other); s
            }
        }
        impl $trait < &$t > for &$t {
            type Output = $t;
            fn $fn (self, other:&$t) -> $t {
                let mut s = *self; s.$f(&other); s
            }
        }
    }
}

macro_rules! binary_assign_op {
    {$t:ty, $trait:ident, $fn:ident, $f:ident} => {
        impl $trait < $t > for $t {
            fn $fn (&mut self, other:$t) {
                self.$f(&other);
            }
        }
        impl $trait < &$t > for $t {
            fn $fn (&mut self, other:&$t) {
                self.$f(other);
            }
        }
    }
}

use std::ops::*;
macro_rules! fp_arith {
    {$t: ty} => {
        binary_op! {$t, Add, add, do_add}
        binary_assign_op! {$t, AddAssign, add_assign, do_add}
        binary_op! {$t, Sub, sub, do_sub}
        binary_assign_op! {$t, SubAssign, sub_assign, do_sub}
        binary_op! {$t, Mul, mul, do_mul}
        binary_assign_op! {$t, MulAssign, mul_assign, do_mul}
        binary_op! {$t, Div, div, do_div}
        binary_assign_op! {$t, DivAssign, div_assign, do_div}
        binary_op! {$t, Rem, rem, do_rem}
        binary_assign_op! {$t, RemAssign, rem_assign, do_rem}
        binary_op! {$t, BitAnd, bitand, do_bit_and}
        binary_assign_op! {$t, BitAndAssign, bitand_assign, do_bit_and}
        binary_op! {$t, BitOr, bitor, do_bit_or}
        binary_assign_op! {$t, BitOrAssign, bitor_assign, do_bit_or}
        binary_op! {$t, BitXor, bitxor, do_bit_xor}
        binary_assign_op! {$t, BitXorAssign, bitxor_assign, do_bit_xor}
    }
}

fp_arith! {FixedPoint_i32_16}
fp_arith! {UnsignedRaw3232}
fp_arith! {SignedRaw3232}
