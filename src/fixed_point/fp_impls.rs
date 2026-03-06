use num_traits::ConstOne;

use super::{FPType, HowIsFixedPoint};

macro_rules! fp_impl{
    {$t:ty, $nb:expr, $nf:expr} => {
        impl HowIsFixedPoint<$t> for FPType<$t, $nf> {
            const NB: usize = $nb;
            const NB_FRAC: usize = $nf;
            const ONE: $t = < $t >::ONE << Self::NB_FRAC;
        }
     }
}

fp_impl! {i8, 8, 1}
fp_impl! {i8, 8, 2}
fp_impl! {i8, 8, 3}
fp_impl! {i8, 8, 4}
fp_impl! {i8, 8, 5}
fp_impl! {i8, 8, 6}

fp_impl! {i32, 32, 1}
fp_impl! {i32, 32, 2}
fp_impl! {i32, 32, 3}
fp_impl! {i32, 32, 4}
fp_impl! {i32, 32, 5}
fp_impl! {i32, 32, 6}
fp_impl! {i32, 32, 7}
fp_impl! {i32, 32, 8}
fp_impl! {i32, 32, 9}
fp_impl! {i32, 32, 10}
fp_impl! {i32, 32, 11}
fp_impl! {i32, 32, 12}
fp_impl! {i32, 32, 13}
fp_impl! {i32, 32, 14}
fp_impl! {i32, 32, 15}
fp_impl! {i32, 32, 16}
fp_impl! {i32, 32, 17}
fp_impl! {i32, 32, 18}
fp_impl! {i32, 32, 19}
fp_impl! {i32, 32, 20}
fp_impl! {i32, 32, 21}
fp_impl! {i32, 32, 22}
fp_impl! {i32, 32, 23}
fp_impl! {i32, 32, 24}
fp_impl! {i32, 32, 25}
fp_impl! {i32, 32, 26}
fp_impl! {i32, 32, 27}
fp_impl! {i32, 32, 28}
fp_impl! {i32, 32, 29}
fp_impl! {i32, 32, 30}

fp_impl! {i64, 64, 60}
