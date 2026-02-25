pub trait Trig: Sized {
    fn sincos_first_quad(self) -> (Self, Self);
}

pub fn apply_rotation_table(
    mut value: i32,
    mut vector: (i32, i32),
    value_table: &[i32],
    rotate_table_pow2: &[u8],
) -> (i32, (i32, i32)) {
    for (v, p) in value_table.iter().zip(rotate_table_pow2.iter()) {
        if value >= 0 {
            value -= *v;
            vector = (vector.0 - (vector.1 >> (*p)), vector.1 + (vector.0 >> (*p)));
        } else {
            value += *v;
            vector = (vector.0 + (vector.1 >> (*p)), vector.1 - (vector.0 >> (*p)));
        }
    }
    (value, vector)
}

pub const NEG_POW2_I32_30: &[u8] = &[
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30,
];
pub const ATAN_ANGLES_I32_30: &[i32] = &[
    0x3243f6a9, 0x1dac6705, 0xfadbafd, 0x7f56ea7, 0x3feab77, 0x1ffd55c, 0xfffaab, 0x7fff55,
    0x3fffeb, 0x1ffffd, 0x100000, 0x80000, 0x40000, 0x20000, 0x10000, 0x8000, 0x4000, 0x2000,
    0x1000, 0x800, 0x400, 0x200, 0x100, 0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1,
];
pub const HALF_PI_U32_30: u32 = 0x6487ed51;
pub const COS_SCALE_U32_30: u32 = 0x26dd3b6a;

impl Trig for i32 {
    fn sincos_first_quad(self) -> (i32, i32) {
        apply_rotation_table(
            self,
            (COS_SCALE_U32_30 as i32, 0),
            ATAN_ANGLES_I32_30,
            NEG_POW2_I32_30,
        )
        .1
    }
}

#[test]
fn test_sincos() {
    fn i32_30(angle: f32) -> i32 {
        (angle * 2.0_f32.powi(30)) as i32
    }
    fn from_i32_30(angle: i32) -> f32 {
        (angle as f32) / 2.0_f32.powi(30)
    }
    for t in crate::utils::float_iter::<f32>(10000) {
        let angle_f32: f32 = std::f32::consts::FRAC_PI_2 * t;
        let (s_f32, c_f32) = angle_f32.sin_cos();
        let angle: i32 = i32_30(angle_f32);
        let (c_i32, s_i32) = angle.sincos_first_quad();
        eprintln!(
            "{angle:08x} {s_i32:08x} {c_i32:08x} {angle_f32} {s_f32} {c_f32} {} {}",
            from_i32_30(s_i32),
            from_i32_30(c_i32),
        );
        assert!(
            (from_i32_30(s_i32) - s_f32).abs() < 1.0E-7,
            "{t} s_i32={s_i32:08x} s_f32={:08x} ({} {s_f32})",
            i32_30(s_f32),
            from_i32_30(s_i32)
        );
        assert!(
            (from_i32_30(c_i32) - c_f32).abs() < 1.0E-7,
            "{t} c_i32={c_i32:08x} c_f32={:08x} ({} {c_f32})",
            i32_30(c_f32),
            from_i32_30(c_i32)
        );
    }
}
