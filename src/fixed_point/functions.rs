/// Takes a value that is u64_0
///
/// Returns a value that is u64_32
///
/// i.e. `x^2 -> x * 2^32`
pub fn sqrt_u64(value: u64) -> u64 {
    let value = (value as u128) << 64;
    let mut est_sqrt: u128 = 0;
    let mut est_sqrt_sq: u128 = 0;
    for i in (0..64).rev() {
        let new_est_sqrt_sq = est_sqrt_sq + (est_sqrt << (i + 1)) + (1 << (2 * i));
        if new_est_sqrt_sq <= value {
            est_sqrt_sq = new_est_sqrt_sq;
            est_sqrt += (1_u64 << i) as u128;
        }
    }
    est_sqrt as u64
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
        let f2 = (v as f64);
        let s = sqrt_u64(v);
        let f = (s as f64) / ((1_u64 << 32) as f64);
        let err = f2.sqrt() - f;
        eprintln!("{v:016x} sqrt -> {:016x} ({f:.6} err:{err})", s);
        assert_eq!(e, s);
    }
}

pub fn apply_rotation_table<const N: usize>(
    mut value: i32,
    mut vector: (i32, i32),
    value_table: &[i32],
    rotate_table_pow2: &[u8],
) -> (i32, (i32, i32)) {
    for (v, p) in value_table.iter().zip(rotate_table_pow2.iter()) {
        for _ in 0..N {
            if value >= 0 {
                value -= *v;
                vector = (vector.0 - (vector.1 >> (*p)), vector.1 + (vector.0 >> (*p)));
            } else {
                value += *v;
                vector = (vector.0 + (vector.1 >> (*p)), vector.1 - (vector.0 >> (*p)));
            }
        }
    }
    (value, vector)
}

pub fn apply_rotation_table_inv<const N: usize>(
    mut vector: (i32, i32),
    atan_angles: &[i32],
    tan_powers: &[u8],
) -> (i32, (i32, i32)) {
    let mut angle = 0;
    for (atan, power) in atan_angles.iter().zip(tan_powers.iter()) {
        for _ in 0..N {
            if vector.1 >= 0 {
                angle += *atan;
                vector = (
                    vector.0 - (vector.1 >> (*power)),
                    vector.1 + (vector.0 >> (*power)),
                );
            } else {
                angle -= *atan;
                vector = (
                    vector.0 + (vector.1 >> (*power)),
                    vector.1 - (vector.0 >> (*power)),
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
        super::apply_rotation_table::<1>(
            angle,
            (COS_SCALE_U32_30 as i32, 0),
            ATAN_ANGLES_I32_30,
            NEG_POW2_I32_30,
        )
        .1
    }

    pub fn atan2(y: i32, x: i32) -> i32 {
        super::apply_rotation_table_inv::<2>((x, y), ATAN_ANGLES_I32_30, NEG_POW2_I32_30).0
    }

    pub fn asin(s: i32) -> i32 {
        let c = 0; // sqrt(1-self^2)
        super::apply_rotation_table_inv::<2>((c, s), ATAN_ANGLES_I32_30, NEG_POW2_I32_30).0
    }

    pub fn acos(c: i32) -> i32 {
        let s = 0; // sqrt(1-self^2)
        super::apply_rotation_table_inv::<2>((c, s), ATAN_ANGLES_I32_30, NEG_POW2_I32_30).0
    }
    fn sqrt(x_sq: i32) -> Option<i32> {
        if x_sq < 0 {
            None
        } else {
            // Use a u64_60 with top bit clear, i.e. max of 2<<60
            let x_sq_64_60 = (x_sq as u64) << 30;
            let x_est_u64_32 = super::sqrt_u64(x_sq_64_60);
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
