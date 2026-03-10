use crate::bignum::UIntN;

use super::{utils, BackingKind, TestKind};

use super::calculate_constants;
use num_traits::FloatConst;

fn hex_array(data: &[u64]) -> String {
    let mut result = format!("[");
    let mut has_nonzero = false;
    for i in data {
        if has_nonzero || *i != 0 {
            result.push_str(&format!("{:#016x}, ", *i));
            has_nonzero = true;
        }
    }
    result.push_str("]");
    result
}

fn calcuate_constants<const N: usize>(frac_bits: u32, display_frac_bits: u32) {
    let one_dbl_frac = UIntN::<N>::with_bit_set(frac_bits * 2);
    let one = UIntN::<N>::with_bit_set(frac_bits);
    let two = one + one;
    let sqrt_two = (two << frac_bits).sqrt();
    let pi = calculate_constants::pi_scaled_by::<N>(frac_bits);
    let e = calculate_constants::e_scaled_by::<N>(frac_bits);
    let frac_2_pi = (one_dbl_frac + one_dbl_frac) / pi;
    let frac_pi_3 = (pi << frac_bits) / (one + one + one);
    let sqrt_pi = calculate_constants::pi_scaled_by::<N>(frac_bits * 2).sqrt();
    let frac_1_sqrt_pi = one_dbl_frac / sqrt_pi;
    let ln_two = calculate_constants::ln_two_scaled_by::<N>(frac_bits);
    let ln_ten = calculate_constants::ln_one_plus_two_neg_power::<N>(frac_bits, 2, None)
        + ln_two
        + ln_two
        + ln_two;
    let log2_e = one_dbl_frac / ln_two;
    let log10_e = one_dbl_frac / ln_ten;
    let log2_ten = (log2_e << frac_bits) / log10_e;
    let log10_two = (log10_e << frac_bits) / log2_e;

    let display_shift = (frac_bits - display_frac_bits) as usize;
    eprintln!(
        "const SQRT_2 = {}",
        &hex_array((sqrt_two >> (display_shift + 1)).raw())
    );
    eprintln!(
        "const LN_2 = {}",
        &hex_array((ln_two >> (display_shift)).raw())
    );
    eprintln!(
        "const LN_10 = {}",
        &hex_array((ln_ten >> (display_shift + 2)).raw())
    );
    eprintln!(
        "const LOG2_E = {}",
        &hex_array((log2_e >> (display_shift + 1)).raw())
    );
    eprintln!(
        "const LOG2_10 = {}",
        &hex_array((log2_ten >> (display_shift + 2)).raw())
    );
    eprintln!(
        "const LOG10_E = {}",
        &hex_array((log10_e >> (display_shift - 1)).raw())
    );
    eprintln!(
        "const LOG10_2 = {}",
        &hex_array((log10_two >> (display_shift - 1)).raw())
    );
    eprintln!(
        "const PI = {}",
        &hex_array((pi >> (display_shift + 2)).raw())
    );
    eprintln!(
        "const FRAC_2_PI = {}",
        &hex_array((frac_2_pi >> display_shift).raw())
    );
    eprintln!(
        "const FRAC_PI_3 = {}",
        &hex_array((frac_pi_3 >> (display_shift + 1)).raw())
    );
    eprintln!(
        "const FRAC_1_SQRT_PI = {}",
        &hex_array((frac_1_sqrt_pi >> (display_shift)).raw())
    );
    eprintln!("const E = {}", &hex_array((e >> (display_shift + 2)).raw()));
}

#[test]
fn calculate_constants_256() {
    // Display of 256 requires 4*u64 + integer bits
    // Double that provides for 1/N
    // But for Sqrt(1/N) we need four times it
    calcuate_constants::<20>(512, 256);
    // assert!(false);
}
#[test]
fn calcuation_of_constants() {
    let pi_126 = calculate_constants::pi_scaled_by::<3>(126);
    eprintln!("{pi_126:#034x} {pi_126}");
    assert_eq!(
        pi_126.integer_decode(),
        (0xc90f_daa2_2168_c234, 63, 1),
        "Error in value of PI"
    );

    let e_126 = calculate_constants::e_scaled_by::<2>(126);
    eprintln!("{e_126:#034x} {e_126}");
    assert_eq!(
        e_126.integer_decode(),
        (0xadf8_5458_a2bb_4a9a, 63, 1),
        "Error in value of e"
    );

    let ln_two_128 = calculate_constants::ln_two_scaled_by::<3>(128);
    eprintln!("{ln_two_128:#034x} {ln_two_128}");
    assert_eq!(
        ln_two_128.integer_decode(),
        (0xb172_17f7_d1cf_79ab, 63, 1),
        "Error in value of ln(2)"
    );

    let ln_five_quarters_128 = calculate_constants::ln_one_plus_two_neg_power::<3>(128, 2, None);
    let ln_ten_128 = ln_five_quarters_128 + ln_two_128 + ln_two_128 + ln_two_128;
    eprintln!("ln(10) = {ln_ten_128:#034x} {ln_ten_128}");
    assert_eq!(
        ln_ten_128.integer_decode(),
        (0x935d_8ddd_aaa8_ac16, 65, 1),
        "Error in value of ln(10)"
    );

    let one_257 = UIntN::<5>::with_bit_set(257);

    let ln_two_128: UIntN<5> = (&ln_two_128).try_into().unwrap();
    let log2_e_129 = one_257 / ln_two_128;
    eprintln!("log2(e) = {log2_e_129:#035x} {log2_e_129}");
    assert_eq!(
        log2_e_129.integer_decode(),
        (0xb8aa_3b29_5c17_f0bb, 65, 1),
        "Error in value of log2(e)"
    );

    let ln_ten_128: UIntN<5> = (&ln_ten_128).try_into().unwrap();
    let log10_e_129 = one_257 / ln_ten_128;
    eprintln!("log10(e) = {log10_e_129:#035x} {log10_e_129}");
    assert_eq!(
        log10_e_129.integer_decode(),
        (0xde5b_d8a9_3728_7195, 63, 1),
        "Error in value of log10(e)"
    );

    let log2_ten_129 = (log2_e_129 << 129_u32) / log10_e_129;
    eprintln!("log2(10) = {log2_ten_129:#035x} {log2_ten_129}");
    assert_eq!(
        log2_ten_129.integer_decode(),
        (0xd49a_784b_cd1b_8afe, 66, 1),
        "Error in value of log2(10)"
    );

    let log10_two_129 = (log10_e_129 << 129_u32) / log2_e_129;
    eprintln!("log10(2) = {log10_two_129:#035x} {log10_two_129}");
    assert_eq!(
        log10_two_129.integer_decode(),
        (0x9a20_9a84_fbcf_f798, 63, 1),
        "Error in value of log10(2)"
    );

    eprintln!("atan table - needs PI/2 as its first entry");
    for i in 1..20 {
        let x = calculate_constants::atan_two_neg_power::<4>(128, i);
        let (m, e, s) = x.integer_decode();
        let x_f = (m as f64) * (s as f64) * (2.0_f64.powi(e as i32 - 127));
        let t_f = ((0.5_f64).powi(i as i32)).atan();
        let diff = (t_f - x_f).abs();
        eprintln!("{i} : {x:#034x} {x_f} {t_f}");
        assert!(diff < 1E-10, "Difference in atan(2^-{i}) was too large");
    }

    // assert!(false, "Force fail");
}

/// An explicit test that requires add to work, that tests 0, 1, 2, 4, 6 and half, using a small amount of raw
///
/// This will not try to test 2, or 4 and 6, if the number of integer bits is too small
pub fn simple_constants<T, I, const N: usize>()
where
    T: TestKind<I, N>,
    I: BackingKind,
{
    let int_bits = utils::int_bits(&T::ZERO);

    let z = T::ZERO;
    let o = T::ONE;

    eprintln!(
        "Simple_constants test on type {} with {N} fractional bits {int_bits} integer bits",
        std::any::type_name::<I>()
    );

    assert!(z.is_zero(), "Zero must be zero");
    assert!(o.is_one(), "One must be one");

    eprintln!("Zero {z:?} One {o:?}");

    assert_ne!(z, o, "One is not zero");

    assert_eq!(z, T::from_u8(0).unwrap(), "Zero must be zero");
    assert_eq!(o, T::from_u8(1).unwrap(), "One must be one");

    assert_eq!(z.integer_decode().0, 0, "Zero encoded as 0");
    assert_eq!(o.integer_decode().2, 1, "One encoded as positive");
    assert_eq!(
        o.integer_decode().0 >> (-o.integer_decode().1),
        1,
        "One encoded as 1<<N"
    );

    if int_bits > 1 {
        let t = o + o;
        assert_eq!(t, T::from_u8(2).unwrap(), "Two must be two");

        eprintln!("Zero {z:?} One {o:?} Two {t:?}");

        assert!(t > o, "2>1");
        assert!(o < t, "1<2");

        if int_bits > 2 {
            let f = t + t;
            assert_eq!(f, T::from_u8(4).unwrap(), "Four must be four");

            assert!(f > t, "4>2");
            assert!(t < f, "2<4");

            let s = t + f;
            assert_eq!(s, T::from_u8(6).unwrap(), "Six must be six");

            assert!(s > f, "6>4");
            assert!(f < s, "4<6");
        }
    }
}

/// An explicit test that requires sub and mul to work fully
///
/// With the arithmetic operations this can test that (e.g.) 1/PI * PI == 1
///
/// This tests TAU, PI, PI/2, PI/3, PI/4, PI/6, PI/8, 1/PI, 2/PI, 2/PI.sqrt()
///
/// SQRT(2), 1/SQRT(2)
///
/// E, LN2, LN10
///
/// LOG2(E), LOG2(10)
///
/// LOG10(E), LOG10(2)
pub fn float_constants<T, I, const N: usize>()
where
    T: TestKind<I, N> + num_traits::FloatConst + num_traits::Float,
    I: BackingKind,
{
    let int_bits = utils::int_bits(&T::ZERO);
    eprintln!(
        "Float_constants test on type {} with {N} fractional bits {int_bits} integer bits",
        std::any::type_name::<I>()
    );

    let mut half = T::ONE;
    half = half / (half + half);
    assert_eq!(half + half, T::ONE);

    let half_pi = T::FRAC_PI_2();

    if int_bits >= 2 {
        let pi = T::PI();
        eprintln!(
            "PI :{:0x} << {}",
            pi.integer_decode().0,
            pi.integer_decode().1
        );

        assert!(
            utils::nearly_equal(pi, half_pi + half_pi, 2),
            "PI - 2*HALF_PI must be different only in rounding"
        );

        eprintln!("PI minus 3 :{:?}", pi - T::from_u8(3).unwrap());
        eprintln!("PI minus 3.5 :{:?}", pi - half - T::from_u8(3).unwrap());
        assert!(
            pi > T::from_u8(3).unwrap(),
            "PI is supposed to be more than 3!"
        );
        assert!(
            pi - half < T::from_u8(3).unwrap(),
            "PI is supposed to be less than 3.5!"
        );

        // The test value is only 48 bits of fraction, so check as many of those as possible
        //
        // This assumes no rounding, which is poor...
        assert!(
            utils::matches_float(pi, f64::PI()),
            "Top bits (up to 48) should match known good value for PI"
        );
    }
    if int_bits >= 3 {
        let tau = T::TAU();
        eprintln!(
            "TAU :{:0x} << {}",
            tau.integer_decode().0,
            tau.integer_decode().1
        );
        assert!(tau > T::from_u8(6).unwrap(), "TAU is more than 6!");
        assert!(tau - half < T::from_u8(6).unwrap(), "TAU is less than 6.5!");

        assert!(
            utils::nearly_equal(tau, half_pi + half_pi + half_pi + half_pi, 4),
            "TAU - 4*HALF_PI must be different only in rounding"
        );
        assert!(
            utils::matches_float(tau, f64::TAU()),
            "Top bits (up to 48) should match known good value for TAU"
        );
    }

    eprintln!(
        "PI/2 :{:0x} << {}",
        half_pi.integer_decode().0,
        half_pi.integer_decode().1
    );
    assert!(half_pi > T::from_u8(1).unwrap(), "PI/2 is more than 1!");
    assert!(
        half_pi - half > T::from_u8(1).unwrap(),
        "PI/2 is more than 1.5!"
    );
    if int_bits >= 2 {
        assert!(half_pi < T::from_u8(2).unwrap(), "PI/2 is less than 2!");
    }
    assert!(
        utils::matches_float(half_pi, f64::FRAC_PI_2()),
        "Top bits (up to 48) should match known good value for half PI"
    );

    if int_bits >= 2 {
        let r_pi = T::FRAC_1_PI();
        let pi = T::PI();
        let r = r_pi * pi;
        dbg!(r);
        assert!(utils::nearly_equal(r, T::ONE, 3), "1/PI * PI should be 1");

        let r = half_pi * T::FRAC_2_PI();
        dbg!(r);
        assert!(utils::nearly_equal(r, T::ONE, 3), "2/PI * PI/2 should be 1");

        assert!(
            utils::matches_float(T::FRAC_2_PI(), f64::FRAC_2_PI()),
            "Top bits (up to 48) should match known good value for 2/PI"
        );
    }

    if int_bits >= 3 {
        let r_pi = T::FRAC_1_PI();
        let r = r_pi * T::TAU() * half;
        assert!(
            utils::nearly_equal(r, T::ONE, 3),
            "1/PI * TAU * half should be 1"
        );
    }

    let r_pi = T::FRAC_1_PI();
    let r = r_pi * half_pi;
    assert!(utils::nearly_equal(r, half, 2), "1/PI * PI/2 should be 1/2");

    assert!(
        utils::nearly_equal(T::FRAC_PI_4(), half * T::FRAC_PI_2(), 2),
        "PI/4 should be PI/2 * 1/2"
    );

    assert!(
        utils::nearly_equal(T::FRAC_PI_8(), half * T::FRAC_PI_4(), 2),
        "PI/8 should be PI/4 * 1/2"
    );

    assert!(
        utils::nearly_equal(T::FRAC_PI_6(), half * T::FRAC_PI_3(), 2),
        "PI/6 should be PI/3 * 1/2"
    );

    assert!(
        utils::nearly_equal(
            T::FRAC_PI_6() + T::FRAC_PI_6() + T::FRAC_PI_6(),
            T::FRAC_PI_2(),
            3
        ),
        "PI/6 + PI/6 + PI/6  should be PI/2"
    );

    eprintln!("{:?}", T::FRAC_2_SQRT_PI() * T::FRAC_2_SQRT_PI() * T::PI());
    assert!(
        utils::nearly_equal(
            T::FRAC_2_SQRT_PI() * half * half * T::FRAC_2_SQRT_PI() * T::PI(),
            T::ONE,
            8
        ),
        "2/PI.sqrt() /4 * 2/PI.sqrt() should be PI"
    );

    assert!(
        utils::nearly_equal(T::SQRT_2() * half, T::FRAC_1_SQRT_2(), 2),
        "sqrt(2)/2 should be 1/sqrt(2)",
    );

    assert!(
        utils::nearly_equal(T::SQRT_2() * T::FRAC_1_SQRT_2(), T::ONE, 2),
        "sqrt(2) * 1/sqrt(2) should be 1",
    );
    assert!(
        utils::nearly_equal(T::FRAC_1_SQRT_2() * T::FRAC_1_SQRT_2(), half, 2),
        "1/sqrt(2)*1/sqrt(2) should be half",
    );

    assert!(
        utils::matches_float(T::E(), f64::E()),
        "Top bits (up to 48) should match known good value for e"
    );
    assert!(
        utils::matches_float(T::LN_2(), f64::LN_2()),
        "Top bits (up to 48) should match known good value for ln(2)"
    );
    assert!(
        utils::matches_float(T::LN_10(), f64::LN_10()),
        "Top bits (up to 48) should match known good value for ln(10)"
    );

    assert!(
        utils::nearly_equal(T::LOG2_E() * T::LN_2(), T::ONE, 3),
        "log2_E"
    );
    assert!(
        utils::nearly_equal(T::LOG2_10() * T::LN_2(), T::LN_10(), 3),
        "log2_10"
    );

    assert!(
        utils::nearly_equal(T::LOG10_E() * T::LN_10(), T::ONE, 3),
        "log10_E"
    );
    assert!(
        utils::nearly_equal(T::LOG10_2() * T::LN_10(), T::LN_2(), 3),
        "log10_2"
    );
}
