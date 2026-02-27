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
    + std::ops::MulAssign<Self>
    + std::ops::Not<Output = Self>
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
        + std::ops::MulAssign<Self>
        + std::ops::Not<Output = Self>
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
pub fn sqrt_dbl<F: UsefulUInt>(value: F::Dbl) -> F {
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

/// Takes a value `x_sq` that is an `N` bit unsigned integer such as u64
///
/// Returns a value that is the same type, whose value is `x * 2^N`
///
/// i.e. `x^2 -> x * 2^32` for u64 (or u64_32)
///
/// If 'x_sq' is a fixed point number with `2B` fractional bits, such as u64_60, then x_sq is really `y_sq*2^(2B)`; x is really `y*2B`, and the return result
/// is `y * 2^(32+B)` (i.e. u64_62)
///
/// For example, a u32_24 input produces a u32_28 output; a u64_0 (pure integer) produces a u64_32 output
#[inline]
pub fn sqrt<F: UsefulUInt>(x_sq: F) -> F {
    sqrt_dbl(x_sq.as_dbl_upper())
}

/// Given a pure fractional value of type `F` (i.e. one is F::MAX+1), return sqrt(1-x^2)
///
/// If x is zero then this will return zero
///
/// This works fine if F is `i64 as u64`
pub fn sqrt_max_minus_x_sq<F: UsefulUInt>(x: F) -> F {
    if x.is_zero() {
        F::ZERO
    } else {
        let mut x_sq = x.as_dbl();
        x_sq *= x_sq;
        eprintln!("x {x:08x} x_sq:{x_sq:016x}");
        sqrt_dbl::<F>((!x_sq) + F::Dbl::ONE)
    }
}

/// Takes a value that is a fixed point value, and returns 1/sqrt()
///
/// Returns a value that is the same type - if input is I.F fixed point then result is (F/2).(I+F/2) fixed point.
///
/// If the input value is considered to be u64_32 then the result is u64_48, as
/// the smallest input value (2^-32) squarerooted (2^-16) has a reciprocal of
/// 2^16, hence the largest value of the result is 2^16 (actually, 2^16-1 in the
/// implementation).
pub fn recip_sqrt<F: UsefulUInt>(x: F) -> F {
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
    let one = F::Dbl::ONE << F::NB;
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
                return est_r_sqrt;
            }
        }
    }
    est_r_sqrt
}

/// Apply a log table to find log(value/one), where 0<value/one<2
///
/// If F is u64_60 then `value`, `one` and the log table (each of 1+2^power) should be in those units.
///
/// The log table can be in any base; normally N is 1,
pub fn apply_log_table<F: UsefulInt, const N: usize>(
    value: F,
    one: F,
    log_table: &[F],
    pow2_table: &[u8],
) -> F {
    let mut log = F::ZERO;
    let mut exp = one;
    for (l, p) in log_table.iter().zip(pow2_table.iter().copied()) {
        for _ in 0..N {
            let new_exp = exp + (exp >> p);
            let cmp = new_exp.cmp(&value);
            if cmp != std::cmp::Ordering::Greater {
                exp = new_exp;
                log += *l;
                if cmp == std::cmp::Ordering::Equal {
                    return log;
                }
            }
        }
    }
    log
}

/// Apply a log table to find -log(one/value), where 0<one/value<2*one
///
pub fn apply_m_log_table<F: UsefulInt, const N: usize>(
    mut value: F,
    one: F,
    log_table: &[F],
    pow2_table: &[u8],
) -> F {
    let mut log = F::ZERO;
    let mut exp = one;
    for (l, p) in log_table.iter().zip(pow2_table.iter().copied()) {
        for _ in 0..N {
            let new_exp = exp - (exp >> (p + 1));
            let cmp = new_exp.cmp(&value);
            if cmp != std::cmp::Ordering::Less {
                exp = new_exp;
                log += *l;
                if cmp == std::cmp::Ordering::Equal {
                    return log;
                }
            }
        }
    }
    log
}

/// Apply a log table to find log(value)
pub fn apply_inv_log_table<F: UsefulInt, const N: usize>(
    mut value: F,
    one: F,
    log_table: &[F],
    pow2_table: &[u8],
) -> F {
    let mut log = F::ZERO;
    let mut exp = one;
    for (l, p) in log_table.iter().zip(pow2_table.iter().copied()) {
        for _ in 0..N {
            let new_log = log + *l;
            let cmp = new_log.cmp(&value);
            if cmp != std::cmp::Ordering::Greater {
                log = new_log;
                exp += exp >> p;
                if cmp == std::cmp::Ordering::Equal {
                    return exp;
                }
            }
        }
    }
    exp
}

pub fn apply_inv_m_log_table<F: UsefulInt, const N: usize>(
    mut value: F,
    one: F,
    log_table: &[F],
    pow2_table: &[u8],
) -> F {
    let mut log = F::ZERO;
    let mut exp = one;
    for (l, p) in log_table.iter().zip(pow2_table.iter().copied()) {
        for _ in 0..N {
            let new_log = log + *l;
            let cmp = new_log.cmp(&value);
            if cmp != std::cmp::Ordering::Greater {
                log = new_log;
                exp -= exp >> (p + 1);
                if cmp == std::cmp::Ordering::Equal {
                    return exp;
                }
            }
        }
    }
    exp
}

pub fn apply_rotation_table<F: UsefulInt, const N: usize>(
    mut angle: F,
    mut vector: (F, F),
    value_table: &[F],
    rotate_table_pow2: &[u8],
) -> (F, (F, F)) {
    // vector.0 should *always* be positive; vector.1 may be negative
    for (v, p) in value_table.iter().zip(rotate_table_pow2.iter().copied()) {
        for _ in 0..N {
            // This loop ought to be able to bail out early
            if angle >= F::ZERO {
                angle -= *v;
                vector = (vector.0 - (vector.1 >> p), vector.1 + (vector.0 >> p));
            } else {
                angle += *v;
                vector = (vector.0 + (vector.1 >> p), vector.1 - (vector.0 >> p));
            }
        }
    }
    (angle, vector)
}

pub fn apply_rotation_table_inv<F: UsefulInt, const N: usize>(
    mut vector: (F, F),
    atan_angles: &[F],
    tan_powers: &[u8],
) -> (F, (F, F)) {
    let mut angle = F::ZERO;
    for (atan, power) in atan_angles.iter().zip(tan_powers.iter().copied()) {
        for _ in 0..N {
            //eprintln!(
            //    "angle {angle:x} atan {atan:x} vector:{:x},{:x}",
            //    vector.0, vector.1,
            //);
            // Note that ths loop cannot bail early, as the angle *must* go
            // through the same number of iterations as it will be scaled
            // dependent on the iterations
            if vector.1 <= F::ZERO {
                angle -= *atan;
                vector = (
                    vector.0 - (vector.1 >> power),
                    vector.1 + (vector.0 >> power),
                );
            } else {
                angle += *atan;
                vector = (
                    vector.0 + (vector.1 >> power),
                    vector.1 - (vector.0 >> power),
                );
            }
        }
    }
    (angle, vector)
}

pub mod i32_28 {
    use crate::fixed_point::functions::{
        apply_inv_log_table, apply_inv_m_log_table, apply_log_table, apply_rotation_table,
        apply_rotation_table_inv,
    };

    pub const NEG_POW2_I32_30: &[u8] = &[
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
        25, 26, 27, 28, 29, 30,
    ];
    pub const ATAN_ANGLES_I32_30: &[i32] = &[
        0x3243f6a9, 0x1dac6705, 0xfadbafd, 0x7f56ea7, 0x3feab77, 0x1ffd55c, 0xfffaab, 0x7fff55,
        0x3fffeb, 0x1ffffd, 0x100000, 0x80000, 0x40000, 0x20000, 0x10000, 0x8000, 0x4000, 0x2000,
        0x1000, 0x800, 0x400, 0x200, 0x100, 0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1,
    ];
    pub const NEG_POW2_LOG_I32_30: &[u8] = &[
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
        25, 26, 27, 28, 29, 30,
    ];
    pub const LOG_POWER_I32_30: &[i32] = &[
        0x2c5c85fe, 0x19f323ed, 0xe47fbe4, 0x789c1dc, 0x3e14618, 0x1f829b1, 0xfe0546, 0x7f80aa,
        0x3fe015, 0x1ff803, 0xffe00, 0x7ff80, 0x3ffe0, 0x1fff8, 0xfffe, 0x8000, 0x4000, 0x2000,
        0x1000, 0x800, 0x400, 0x200, 0x100, 0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1,
    ];
    pub const INV_LOG_POWER_I32_30: &[i32] = &[
        0x2c5c85fe, 0x12696211, 0x88bc741, 0x421662d, 0x2082bb1, 0x1020566, 0x8080ac, 0x402015,
        0x200803, 0x100200, 0x80080, 0x40020, 0x20008, 0x10002, 0x8001, 0x4000, 0x2000, 0x1000,
        0x800, 0x400, 0x200, 0x100, 0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1, 0x1,
    ];
    pub const HALF_PI_U32_30: u32 = 0x6487ed51;
    pub const COS_SCALE_U32_30: u32 = 0x26dd3b6a;

    pub const NEG_POW2_LOG_I32_28: &[u8] = &[
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
        25, 26, 27, 28,
    ];
    // log_table[i] = ln(1 + 2^-i),
    pub const LOG_POWER_I32_28: &[i32] = &[
        0xb17217f, 0x67cc8fb, 0x391fef9, 0x1e27077, 0xf85186, 0x7e0a6c, 0x3f8151, 0x1fe02a,
        0xff805, 0x7fe01, 0x3ff80, 0x1ffe0, 0xfff8, 0x7ffe, 0x4000, 0x2000, 0x1000, 0x800, 0x400,
        0x200, 0x100, 0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1, 0x0, 0x0,
    ];
    // inv_log_table[i] = ln(1 - 1/(2^(i + 1))
    pub const INV_LOG_POWER_I32_28: &[i32] = &[
        0xb17217f, 0x49a5884, 0x222f1d0, 0x108598b, 0x820aec, 0x408159, 0x20202b, 0x100805,
        0x80201, 0x40080, 0x20020, 0x10008, 0x8002, 0x4001, 0x2000, 0x1000, 0x800, 0x400, 0x200,
        0x100, 0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1, 0x1, 0x0, 0x0,
    ];
    pub const NEG_POW2_I32_28: &[u8] = &[
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
        25, 26, 27, 28,
    ];
    pub const ATAN_ANGLES_I32_28: &[i32] = &[
        0xc90fdaa, 0x76b19c1, 0x3eb6ebf, 0x1fd5baa, 0xffaade, 0x7ff557, 0x3ffeab, 0x1fffd5,
        0xffffb, 0x7ffff, 0x40000, 0x20000, 0x10000, 0x8000, 0x4000, 0x2000, 0x1000, 0x800, 0x400,
        0x200, 0x100, 0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1,
    ];
    pub const HALF_PI_U32_28: u32 = 0x6487ed51 >> 2;
    pub const COS_SCALE_U32_28: u32 = 0x9b74edb;

    /// Calculate sin and cos of an angle in the range -PI/2 to +PI/2
    ///
    /// Should clamp in range
    pub fn sincos_first_quad(angle: i32) -> (i32, i32) {
        let (c, s) = super::apply_rotation_table::<_, 1>(
            angle,
            (COS_SCALE_U32_28 as i32, 0),
            ATAN_ANGLES_I32_28,
            NEG_POW2_I32_28,
        )
        .1;
        (s, c)
    }

    pub fn atan2(y: i32, x: i32) -> i32 {
        if x > y {
            apply_rotation_table_inv::<_, 2>((x, y), ATAN_ANGLES_I32_28, NEG_POW2_I32_28).0
        } else {
            HALF_PI_U32_28 as i32
                - apply_rotation_table_inv::<_, 2>((y, x), ATAN_ANGLES_I32_28, NEG_POW2_I32_28).0
        }
    }

    pub fn asin(s: i32) -> i32 {
        if s == 0 {
            0
        } else {
            let c = (super::sqrt_max_minus_x_sq((s << 4) as u32) >> 4) as i32;
            atan2(s, c)
        }
    }

    pub fn acos(c: i32) -> i32 {
        if c == 0 {
            HALF_PI_U32_28 as i32
        } else {
            let s = (super::sqrt_max_minus_x_sq((c << 4) as u32) as i32) >> 4;
            atan2(s, c)
        }
    }

    pub fn sqrt(x_sq: i32) -> Option<i32> {
        if x_sq < 0 {
            None
        } else {
            // This is overkill probably - could use just sqrt(x_sq) (which returns a u32_30) and >>2
            //
            // super::sqrt takes a u64_60 and returns a u64_62.
            let x_sq_u64_60 = (x_sq as u64) << 32;
            let x_est_u64_62 = super::sqrt(x_sq_u64_60);
            Some((x_est_u64_62 >> 34) as i32)
        }
    }
    pub fn ln(x: i32) -> Option<i32> {
        Some(super::apply_log_table::<_, 1>(
            x,
            1 << 28,
            LOG_POWER_I32_28,
            NEG_POW2_LOG_I32_28,
        ))
    }

    pub fn ln_recip(x: i32) -> Option<i32> {
        Some(super::apply_m_log_table::<_, 2>(
            x,
            1 << 28,
            INV_LOG_POWER_I32_28,
            NEG_POW2_LOG_I32_28,
        ))
    }

    /// Range supported: ???
    pub fn exp(x: i32) -> Option<i32> {
        Some(super::apply_inv_log_table::<_, 1>(
            x,
            1 << 28,
            LOG_POWER_I32_28,
            NEG_POW2_LOG_I32_28,
        ))
    }

    /// Range supported: ???
    pub fn exp_m(x: i32) -> Option<i32> {
        Some(super::apply_inv_m_log_table::<_, 2>(
            x,
            1 << 28,
            INV_LOG_POWER_I32_28,
            NEG_POW2_LOG_I32_28,
        ))
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

#[test]
fn test_u64_recip() {
    for (v, e) in [
        (2 << 40, 0x00000b504f333f9d),
        (2 << 32, 0x00000b504f333f9de),
        (1 << 32, 1 << 48),
        (1 << 36, 1 << 46),
        (1 << 8, 1 << 60),
        (1 << 0, u64::MAX),
        (2 << 0, 0xb504f333f9de6487),
    ] {
        let f2 = (v as f64) / ((1_u64 << 32) as f64);
        let r_s = recip_sqrt::<u64>(v);
        let r_f = (r_s as f64) / ((1_u64 << 48) as f64);
        let err = r_f - (1.0 / f2.sqrt());
        eprintln!("{v:016x} sqrt -> {r_s:016x} ({r_f:.6} err:{err})");
        assert_eq!(e, r_s);
    }
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
    for (v, e, somv2) in [
        (0, 0x0, 0x0),
        (1, 1 << 32, 0xffffffffffffffff),
        (2, 0x16a09e667, 0xffffffffffffffff),
        (3, 0x1bb67ae85, 0xffffffffffffffff),
        (4, 1 << 33, 0xffffffffffffffff),
        (1 << 16, 1 << 40, 0xffffffffffffffff),
        (2 << 16, 0x16a09e667f3, 0xffffffffffffffff),
        (4 << 16, 2 << 40, 0xffffffffffffffff),
        (1 << 32, 1 << 48, 0xffffffffffffffff),
        (2 << 32, 0x16a09e667f3bc, 0xfffffffffffffffd),
        (4 << 32, 2 << 48, 0xfffffffffffffff7),
        (1 << 56, 1 << 60, 0xffff7fffdfffefff),
        (2 << 56, 0x16a09e667f3bcc90, 0xfffdfffdfffbfff5),
        (1 << 60, 1 << 62, 0xff7fdfeff5f8fabb),
        (2 << 60, 0x5a827999fcef3242, 0xfdfdfbf5e3aaf49a),
        (4 << 60, 0x8000000000000000, 0xf7def58a7a76cd8b),
        (8 << 60, 0xb504f333f9de6484, 0xddb3d742c265539d),
        (0xb504f333f9de6484, 0xd744fccad69d6af4, 0xb504f333f9de6484), // 1/sqrt(2) - so sqrt(1-x^2) = sqrt(1-1/2) = 1/sqrt(2)
        (12 << 60, 0xddb3d742c265539d, 0xa953fd4e97c74dbc),           // 3/4, so sqrt(3)/2
        (0xddb3d742c265539d, 0xee3c1ebb579acafb, 0x8000000000000000), // sqrt(3)/2 - so sqrt(1-x^2) = sqrt(1-3/4) = sqrt(1/4) = 1/2
        (u64::MAX - 1, 0xfffffffffffffffe, 0x00000001ffffffff),
        (u64::MAX, 0xffffffffffffffff, 0x000000016a09e667),
    ] {
        let f2 = v as f64;
        let s = sqrt(v);
        let c = sqrt_max_minus_x_sq(v);
        let f = (s as f64) / ((1_u64 << 32) as f64);
        let err = f2.sqrt() - f;
        eprintln!("{v:016x} sqrt -> {s:016x} ({f:.6} err:{err}) {c:016x}");
        assert_eq!(
            e, s,
            "Mismatch in expected value of square-root for {v:08x}"
        );
        assert_eq!(
            somv2, c,
            "Mismatch in expected value of square-root of one minus x squared for {v:08x}"
        );
    }
}

#[test]
fn test_i32_28_exp() {
    fn i32_28(f: f64) -> i32 {
        (f * 2.0_f64.powi(28)) as i32
    }
    fn from_i32_28(i: i32) -> f64 {
        (i as f64) / 2.0_f64.powi(28)
    }
    let to_float = from_i32_28;
    let from_float = i32_28;

    // for l_f in [1.2, 1.0, 0.75, 0.5, 0.25, 0.125, 0.0625] {
    for i in 1..=500 {
        let l_f = ((i as f64) / 500.0) * 1.4;
        let l = from_float(l_f);
        let e = i32_28::exp(l).expect("Should have managed exp");
        let e_f = to_float(e);

        let e_r = i32_28::exp_m(l).expect("Should have managed exp_m");
        let e_r_f = to_float(e_r);

        if false {
            eprintln!(
                "l:{l_f} e^l:{e_f} {} 1/e^l:{e_r_f} {}",
                l_f.exp(),
                e_r_f * e_f
            );
        }

        assert!(
            (l_f.exp() - e_f).abs() < 1E-5,
            "Exp(l) versus float version should be equal",
        );
        assert!(
            (e_r_f * e_f - 1.0).abs() < 1E-5,
            "Exp(l) * exp(-l) should be 1.0"
        );
    }

    for i in 1..=500 {
        let e_f = ((i as f64) / 500.0) * 3.5 + 1.0;
        //    for e_f in [
        //        1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0,
        //        3.2, 3.4, 3.8, 4.5,
        //    ] {
        let e = from_float(e_f);
        let e_r = from_float(1.0 / e_f);
        let l = i32_28::ln(e).expect("Should have managed ln");
        let l_f = to_float(l);

        let l_r = i32_28::ln_recip(e_r).expect("Should have managed ln_m");
        let l_r_f = to_float(l_r);

        if (false) {
            eprintln!(
                "e:{e_f} ln(e):{l_f} {l:08x} {} {l_r_f} {}",
                e_f.ln(),
                l_r_f - l_f
            );
        }
        assert!(
            (l_f - e_f.ln()).abs() < 1E-5,
            "Ln(e) versus float version should be equal",
        );
        assert!((l_r_f - l_f).abs() < 1E-5, "Ln(e) + ln(1/e) should be 0.0");
    }
    // assert!(false, "Force failure");
}

#[test]
fn test_i32_28_trig() {
    fn i32_28(angle: f64) -> i32 {
        (angle * 2.0_f64.powi(28)) as i32
    }
    fn from_i32_28(angle: i32) -> f64 {
        (angle as f64) / 2.0_f64.powi(28)
    }
    let to_float = from_i32_28;
    let from_float = i32_28;

    /// Vector is overflowing when it is 0.707 * 1<<30 = 0x2d413ccc and we add 0x3243f6a9 (for +atan(0.5)=26.56 degrees twice)
    let x = i32_28::atan2(1 << 27, 1 << 27);
    let x_f = to_float(x);
    eprintln!("{x:x} {x_f} {}", x_f.to_degrees());
    let e = (i32_28::HALF_PI_U32_28 as i32 / 2);
    let err = (e - x).abs();
    assert!(
        err < 4,
        "atan2 of 1,1 is PI/4 (got {x:08x} expected {e:08x}"
    );

    let x = i32_28::atan2(1 << 28, 0);
    let x_f = to_float(x);

    eprintln!("{x:x} {x_f} {}", x_f.to_degrees());
    let e = (i32_28::HALF_PI_U32_28 as i32);
    let err = (e - x).abs();
    assert!(
        err < 4,
        "atan2 of 1,0 is PI/2 (got {x:08x} expected {e:08x}"
    );

    let x = i32_28::atan2(1 << 28, 0x376cf5d0 >> 1);
    let x_f = to_float(x);
    eprintln!("{x:x} {x_f} {}", x_f.to_degrees());
    let e = (i32_28::HALF_PI_U32_28 as i32 / 3);
    let err = (e - x).abs();
    assert!(
        err < 4,
        "atan2 of 1, sqrt(3) is PI/6 (got {x:08x} expected {e:08x}"
    );

    let x = i32_28::asin(0x2d413ccc >> 2);
    let x_f = to_float(x);
    eprintln!("{x:x} {x_f} {}", x_f.to_degrees());
    let e = (i32_28::HALF_PI_U32_28 as i32 / 2);
    let err = (e - x).abs();
    assert!(
        err < 4,
        "asin of 1/sqrt(2) is PI/4 (got {x:08x} expected {e:08x}"
    );

    for i in (-20..=20) {
        let i_f = (i as f64) / 20.0;
        let angle_f = std::f64::consts::FRAC_PI_2 * i_f;
        let angle = from_float(angle_f);
        let (s, c) = i32_28::sincos_first_quad(angle);
        let s_f = to_float(s);
        let c_f = to_float(c);
        let e_s_f = angle_f.sin();
        let e_c_f = angle_f.cos();
        eprintln!(
            "{:8.3} sin/cos {s:x}/{c:x} = {s_f}/{c_f} expected {e_s_f}/{e_c_f}",
            angle_f.to_degrees()
        );
        assert!((e_s_f - s_f).abs() < 1E-4, "Error in sin too large");
        assert!((e_c_f - c_f).abs() < 1E-4, "Error in cos too large");
    }
}
