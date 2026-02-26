use num_traits::{ConstOne, ConstZero, Num};
pub trait Int:
    Copy
    + PartialOrd
    + Ord
    + std::fmt::Display
    + std::fmt::LowerHex
    + Num
    + ConstOne
    + ConstZero
    + std::ops::AddAssign<Self>
    + std::ops::SubAssign<Self>
    + std::ops::Shl<usize, Output = Self>
    + std::ops::Shr<usize, Output = Self>
    + std::ops::Shr<u8, Output = Self>
{
}
impl<T> Int for T where
    T: Copy
        + PartialOrd
        + Ord
        + std::fmt::Display
        + std::fmt::LowerHex
        + Num
        + ConstOne
        + ConstZero
        + std::ops::AddAssign<Self>
        + std::ops::SubAssign<Self>
        + std::ops::Shl<usize, Output = Self>
        + std::ops::Shr<usize, Output = Self>
        + std::ops::Shr<u8, Output = Self>
{
}

pub trait UsefulInt: Int {
    type Unsigned: UsefulUInt;
    type Dbl: Int;
    const NB: usize;
    fn as_dbl(self) -> Self::Dbl;
    fn as_dbl_upper(self) -> Self::Dbl {
        self.as_dbl() << Self::NB
    }
    fn unsigned(self) -> Self::Unsigned;
    fn of_unsigned(v: Self::Unsigned) -> Self;
}

pub trait UsefulUInt: Int {
    type Dbl: Int;
    const NB: usize;
    fn as_dbl(self) -> Self::Dbl;
    fn as_dbl_upper(self) -> Self::Dbl {
        self.as_dbl() << Self::NB
    }
}

impl UsefulUInt for u16 {
    type Dbl = u32;
    const NB: usize = 16;
    fn as_dbl(self) -> u32 {
        self as u32
    }
}

impl UsefulInt for i16 {
    type Unsigned = u16;
    type Dbl = i32;
    const NB: usize = 16;
    fn as_dbl(self) -> i32 {
        self as i32
    }
    fn unsigned(self) -> Self::Unsigned {
        self as u16
    }
    fn of_unsigned(v: Self::Unsigned) -> Self {
        v as Self
    }
}

impl UsefulUInt for u32 {
    type Dbl = u64;
    const NB: usize = 32;
    fn as_dbl(self) -> u64 {
        self as u64
    }
}

impl UsefulInt for i32 {
    type Unsigned = u32;
    type Dbl = i64;
    const NB: usize = 32;
    fn as_dbl(self) -> i64 {
        self as i64
    }
    fn unsigned(self) -> Self::Unsigned {
        self as u32
    }
    fn of_unsigned(v: Self::Unsigned) -> Self {
        v as Self
    }
}

impl UsefulUInt for u64 {
    type Dbl = u128;
    const NB: usize = 64;
    fn as_dbl(self) -> u128 {
        self as u128
    }
}

impl UsefulInt for i64 {
    type Unsigned = u64;
    type Dbl = i128;
    const NB: usize = 64;
    fn as_dbl(self) -> i128 {
        self as i128
    }
    fn unsigned(self) -> Self::Unsigned {
        self as u64
    }
    fn of_unsigned(v: Self::Unsigned) -> Self {
        v as Self
    }
}

/// Takes a value that is an integer
///
/// Returns a value that is the same type but top half integer, bottom half
/// fractional (e.g. u64_32)
///
/// i.e. `x^2 -> x * 2^32` for u64
pub fn sqrt<F: UsefulUInt>(value: F) -> F {
    // FIXME change to as_dbl_upper()
    let value = value.as_dbl() << F::NB;
    let mut est_sqrt = F::ZERO;
    let mut est_sqrt_sq = F::Dbl::ZERO;
    for i in (0..F::NB).rev() {
        let new_est_sqrt_sq =
            est_sqrt_sq + (est_sqrt.as_dbl() << (i + 1)) + (F::Dbl::ONE << (2 * i));
        if new_est_sqrt_sq <= value {
            est_sqrt_sq = new_est_sqrt_sq;
            est_sqrt += (F::ONE << i);
        }
    }
    est_sqrt
}

/// Takes a value that is a fixed point value with half integer half fraction bits (e.g. u64_32)
///
/// The square-root of an I.F fixed point value has I/2 integer bits; the
/// reciprocal could be *any* number of integer and fractional bits, but really
/// it is not worth having many more than that I/2 as the new fractional size, and similarly F/2 as the new integer size
///
/// Howeever, a u64_60 *could* become a u32_30 (as a square root) and then a u32_2 as a reciprocal
///
/// Returns a value that is the same type
///
/// i.e. `x^2 -> x * 2^32` for u64
pub fn recip_sqrt<F: UsefulUInt, const FRAC_BITS: usize>(x: F) -> F {
    // est_r_sqrt = 1 / sqrt(x) => est_r_sqrt^2 * x = 1
    //
    // If est_r_sqrt has '1<<i' added to it, then
    // est_r_sqrt_sq += (est_r_sq<<(i+1)) + (1<<(2*i))
    // x_times_est_r_sqrt_sq += ((x*est_r_sq)<<(i+1)) + (x<<(2*i))
    //
    // So this needs to maintain x*est_r_sq too, and when est_r_sqrt has `1<<i` added to it:
    //
    // x_times_est_r_sq += x<<i
    //
    // This operation is peformed for each bit i in the possible answer, from top bit down, provided est_r_sqrt^2*x remains less than or equal to 1
    let mut est_r_sqrt = F::ZERO; // x is u32_8 => u32_24
    let mut x_times_est_r_sqrt = F::Dbl::ZERO; // x is u32_8 => u64_32
    let mut x_times_est_r_sqrt_sq = F::Dbl::ZERO; // in theory x is u32_8 => u96_56; in practice u64_24
    let x_shift_fbits = x.as_dbl(); // _upper();
    let one = F::Dbl::ONE << (F::NB - FRAC_BITS);
    for i in (0..F::NB).rev() {
        let new_xtesq = {
            if 2 * i < F::NB {
                x_times_est_r_sqrt_sq
                    + (x_times_est_r_sqrt >> (F::NB - i - 1))
                    + (x_shift_fbits >> (F::NB - 2 * i))
            } else {
                x_times_est_r_sqrt_sq
                    + (x_times_est_r_sqrt >> (F::NB - i - 1))
                    + (x_shift_fbits << (2 * i - F::NB))
            }
        };
        let cmp = new_xtesq.cmp(&one);
        if cmp != std::cmp::Ordering::Greater {
            x_times_est_r_sqrt_sq = new_xtesq;
            x_times_est_r_sqrt += x.as_dbl_upper() >> (F::NB - i);
            est_r_sqrt += F::ONE << i;
            if cmp == std::cmp::Ordering::Equal {
                break;
            }
        }
        // eprintln!(
        //     "{i:2} {one:16x} {est_r_sqrt:16x} {x_times_est_r_sqrt:16x} {x_times_est_r_sqrt_sq:16x}",
        // );
    }
    est_r_sqrt
}

#[test]
fn test_u64_recip() {
    for (v, e) in [
        (2 << 40, 0xb504f33),
        (2 << 32, 0xb504f335),
        (1 << 32, 1 << 32),
        (1 << 36, 1 << 30),
        (1 << 8, 1 << 44),
        (1 << 0, 1 << 48),
        (2 << 0, 0x0000b504f3358000),
    ] {
        let f2 = (v as f64) / ((1_u64 << 32) as f64);
        let r_s = recip_sqrt::<u64, 32>(v);
        let r_f = (r_s as f64) / ((1_u64 << 32) as f64);
        let err = r_f - (1.0 / f2.sqrt());
        eprintln!("{v:016x} sqrt -> {r_s:016x} ({r_f:.6} err:{err})");
        assert_eq!(e, r_s);
    }
    assert!(false);
}

#[test]
fn test_u16() {
    for i in 0..16 {
        let v: u16 = 1 << i;
        let f2 = v as f64;
        let s = sqrt(v);
        let f = (s as f64) / ((1_u64 << 8) as f64);
        let err = f2.sqrt() - f;
        eprintln!("{v:016x} sqrt -> {:016x} ({f:.6} err:{err})", s);
        assert!(
            err < 1.0 / 256.0,
            "Err was {err}, max expected was {}",
            1.0 / 256.0
        );
    }
}

#[test]
fn test_u32() {
    for i in 0..32 {
        let v: u32 = 1 << i;
        let f2 = v as f64;
        let s = sqrt(v);
        let f = (s as f64) / ((1_u64 << 16) as f64);
        let err = f2.sqrt() - f;
        eprintln!("{v:016x} sqrt -> {:016x} ({f:.6} err:{err})", s);
        assert!(
            err < 1.0 / 65536.0,
            "Err was {err}, max expected was {}",
            1.0 / 65536.0
        );
    }
}

#[test]
fn test_u64() {
    for (v, e) in [
        (0, 0x0),
        (1, 1 << 32),
        (2, 0x16a09e667),
        (3, 0x1bb67ae85),
        (4, 1 << 33),
        (1 << 16, 1 << 40),
        (2 << 16, 0x16a09e667f3),
        (4 << 16, 2 << 40),
        (1 << 32, 1 << 48),
        (2 << 32, 0x16a09e667f3bc),
        (4 << 32, 2 << 48),
        (1 << 56, 1 << 60),
        (2 << 56, 0x16a09e667f3bcc90),
        (1 << 60, 1 << 62),
        (2 << 60, 0x5a827999fcef3242),
        (u64::MAX - 1, 0xfffffffffffffffe),
        (u64::MAX, 0xffffffffffffffff),
    ] {
        let f2 = v as f64;
        let s = sqrt(v);
        let f = (s as f64) / ((1_u64 << 32) as f64);
        let err = f2.sqrt() - f;
        eprintln!("{v:016x} sqrt -> {:016x} ({f:.6} err:{err})", s);
        assert_eq!(e, s);
    }
}

pub fn apply_rotation_table<F: UsefulInt, const N: usize>(
    mut value: F,
    mut vector: (F, F),
    value_table: &[F],
    rotate_table_pow2: &[u8],
) -> (F, (F, F)) {
    for (v, p) in value_table.iter().zip(rotate_table_pow2.iter().copied()) {
        for _ in 0..N {
            if value >= F::ZERO {
                value -= *v;
                vector = (vector.0 - (vector.1 >> p), vector.1 + (vector.0 >> p));
            } else {
                value += *v;
                vector = (vector.0 + (vector.1 >> p), vector.1 - (vector.0 >> p));
            }
        }
    }
    (value, vector)
}

pub fn apply_rotation_table_inv<F: UsefulInt, const N: usize>(
    mut vector: (F, F),
    atan_angles: &[F],
    tan_powers: &[u8],
) -> (F, (F, F)) {
    let mut angle = F::ZERO;
    for (atan, power) in atan_angles.iter().zip(tan_powers.iter().copied()) {
        for _ in 0..N {
            if vector.1 >= F::ZERO {
                angle += *atan;
                vector = (
                    vector.0 - (vector.1 >> power),
                    vector.1 + (vector.0 >> power),
                );
            } else {
                angle -= *atan;
                vector = (
                    vector.0 + (vector.1 >> power),
                    vector.1 - (vector.0 >> power),
                );
            }
        }
    }
    (angle, vector)
}

pub mod i32_30 {
    pub const NEG_POW2_I32_30: &[u8] = &[
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
        25, 26, 27, 28, 29, 30,
    ];
    pub const ATAN_ANGLES_I32_30: &[i32] = &[
        0x3243f6a9, 0x1dac6705, 0xfadbafd, 0x7f56ea7, 0x3feab77, 0x1ffd55c, 0xfffaab, 0x7fff55,
        0x3fffeb, 0x1ffffd, 0x100000, 0x80000, 0x40000, 0x20000, 0x10000, 0x8000, 0x4000, 0x2000,
        0x1000, 0x800, 0x400, 0x200, 0x100, 0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1,
    ];
    pub const HALF_PI_U32_30: u32 = 0x6487ed51;
    pub const COS_SCALE_U32_30: u32 = 0x26dd3b6a;

    pub fn sincos_first_quad(angle: i32) -> (i32, i32) {
        super::apply_rotation_table::<_, 1>(
            angle,
            (COS_SCALE_U32_30 as i32, 0),
            ATAN_ANGLES_I32_30,
            NEG_POW2_I32_30,
        )
        .1
    }

    pub fn atan2(y: i32, x: i32) -> i32 {
        super::apply_rotation_table_inv::<_, 2>((x, y), ATAN_ANGLES_I32_30, NEG_POW2_I32_30).0
    }

    pub fn asin(s: i32) -> i32 {
        let c = 0; // sqrt(1-self^2)
        super::apply_rotation_table_inv::<_, 2>((c, s), ATAN_ANGLES_I32_30, NEG_POW2_I32_30).0
    }

    pub fn acos(c: i32) -> i32 {
        let s = 0; // sqrt(1-self^2)
        super::apply_rotation_table_inv::<_, 2>((c, s), ATAN_ANGLES_I32_30, NEG_POW2_I32_30).0
    }
    fn sqrt(x_sq: i32) -> Option<i32> {
        if x_sq < 0 {
            None
        } else {
            // Use a u64_60 with top bit clear, i.e. max of 2<<60
            let x_sq_64_60 = (x_sq as u64) << 30;
            let x_est_u64_32 = super::sqrt(x_sq_64_60);
            // x_est_u64_32 is sqrt * 2^31, max of 1.4<<30
            Some((x_est_u64_32 >> 30) as i32)
        }
    }
    /*    fn cbrt(self) -> Self {
            0
        }
        fn recip(self) -> Self {
            0
        }
        fn recip_sqrt(self) -> Self {
            0
        }
        fn hypot(self) -> Self {
            0
        }
    */
}
