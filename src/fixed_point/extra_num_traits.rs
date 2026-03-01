pub trait OverflowingMul: Sized {
    fn overflowing_mul(self, rhs: Self) -> (Self, bool);
}
pub trait CarryingAdd: Sized {
    /// Calculate the subtraction of rhs from self, and a further subtraction of the borrow if set, returning the wrapped value and a boolean indication of overflow
    ///
    /// Since the bool in the result indicates *overflow* not carry out, it must not be used as a carry in to other operations
    fn carrying_add(self, rhs: Self, carry: bool) -> (Self, bool);
}
pub trait BorrowingSub: Sized {
    /// Calculate the sum of self and rhs and a carry in, returning the wrapped value and a boolean indication of overflow
    ///
    /// Since the bool in the result indicates *overflow* not carry out, it must not be used as a carry in to other operations
    fn borrowing_sub(self, rhs: Self, carry: bool) -> (Self, bool);
}

macro_rules! make_signed_int_deps {
    {$t:ty} => {
    impl OverflowingMul for $t {
        fn overflowing_mul(self, rhs: Self) -> (Self, bool) {
            self.overflowing_mul(rhs)
        }
    }
    impl CarryingAdd for $t {
        fn carrying_add(self, rhs: Self, carry: bool) -> (Self, bool) {
            let (r, o) = self.overflowing_add(rhs);
            if !carry {
                (r, o)
            } else {
                let (r2, o2) = r.overflowing_add(1);
                (r2, o | o2)
            }
        }
    }
    impl BorrowingSub for $t {
        fn borrowing_sub(self, rhs: Self, borrow: bool) -> (Self, bool) {
            let (r, o) = self.overflowing_sub(rhs);
            if !borrow {
                (r, o)
            } else {
                let (r2, o2) = r.overflowing_sub(1);
                (r2, o | o2)
            }
        }
    }
    }
}

make_signed_int_deps!(i8);
make_signed_int_deps!(i16);
make_signed_int_deps!(i32);
make_signed_int_deps!(i64);
make_signed_int_deps!(isize);
make_signed_int_deps!(i128);
