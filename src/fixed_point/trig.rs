use crate::fixed_point::functions;

use super::functions::sqrt_u64;

pub trait Roots: Sized {
    fn sqrt(self) -> Option<Self>;
    fn cbrt(self) -> Self;
    fn recip(self) -> Self;
    fn recip_sqrt(self) -> Self;
    fn hypot(self) -> Self;
}

pub trait Trig: Sized {
    fn sincos_first_quad(self) -> (Self, Self);
    fn atan2(self, x: Self) -> Self;
    fn asin(self) -> Self;
    fn acos(self) -> Self;
}

// i32 as i32_30 for now
impl Roots for i32 {
    fn sqrt(self) -> Option<Self> {
        if self < 0 {
            None
        } else {
            // value is u64_60 with top bit clear, i.e. max of 2<<60
            let value = (self as u64) << 30;
            // Treats value as u64, i.e. value * 2^62
            let sqrt_est_u64_32 = sqrt_u64(value);
            // sqrt_est_u64_32 is sqrt * 2^31, max of 1.4<<30
            Some((sqrt_est_u64_32 >> 30) as i32)
        }
    }
    fn cbrt(self) -> Self {
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
}

// i32 as i32_30 for now
impl Trig for i32 {
    fn sincos_first_quad(self) -> (i32, i32) {
        functions::i32_30::sincos_first_quad(self)
    }

    fn atan2(self, x: i32) -> i32 {
        functions::i32_30::atan2(self, x)
    }

    fn asin(self) -> i32 {
        functions::i32_30::asin(self)
    }

    fn acos(self) -> i32 {
        functions::i32_30::acos(self)
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
