use super::{utils, BackingKind, TestKind};

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

    assert_eq!(utils::raw_u64(&z), 0, "Zero encoded as 0");
    assert_eq!(utils::raw_u64(&o), 1 << N, "One encoded as 1<<N");

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

    let mut half = T::ZERO;
    utils::set_raw(&mut half, 1 << (N - 1));
    eprintln!("Half {half:?}");
    assert_eq!(half + half, o, "Two halves added is one");
    assert_eq!(
        utils::raw_u64(&half),
        1 << (N - 1),
        "Half encoded as 1<<(N-1)"
    );

    let mut smallest = T::ZERO;
    utils::set_raw(&mut smallest, 1);
    assert_eq!(utils::raw_u64(&smallest), 1, "Smallest value must be 1");
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
    let mut half = T::ZERO;
    utils::set_raw(&mut half, 1 << (N - 1));
    let pi_raw_u64_48: u64 = 0x3243f6a8885a3;
    let e_raw_u64_48: u64 = 0x2b7e151628aed;
    let ln2_raw_u64_48: u64 = 0xb17217f7d1cf;
    let ln10_raw_u64_48: u64 = 0x24d763776aaa3;

    let half_pi = T::FRAC_PI_2();

    if int_bits >= 2 {
        let pi = T::PI();
        eprintln!("PI :{:0x}", utils::raw_u64(&pi));
        assert!(pi > T::from_u8(3).unwrap(), "PI is more than 3!");
        assert!(pi - half < T::from_u8(3).unwrap(), "PI is less than 3.5!");

        // The test value is only 48 bits of fraction, so check as many of those as possible
        //
        // This assumes no rounding, which is poor...
        if N < 48 {
            assert_eq!(utils::raw_u64(&pi), pi_raw_u64_48 >> (48 - N));
        } else {
            assert_eq!(utils::raw_u64(&pi) >> (N - 48), pi_raw_u64_48);
        }
        assert!(
            utils::nearly_equal(pi, half_pi + half_pi, 2),
            "PI - 2*HALF_PI must be different only in rounding"
        );
    }
    if int_bits >= 3 {
        let tau = T::TAU();
        eprintln!("TAU :{:0x}", utils::raw_u64(&tau));
        assert!(tau > T::from_u8(6).unwrap(), "TAU is more than 6!");
        assert!(tau - half < T::from_u8(6).unwrap(), "TAU is less than 6.5!");

        assert!(
            utils::nearly_equal(tau, half_pi + half_pi + half_pi + half_pi, 4),
            "TAU - 4*HALF_PI must be different only in rounding"
        );
    }

    eprintln!("PI/2 :{:0x}", utils::raw_u64(&half_pi));
    assert!(half_pi > T::from_u8(1).unwrap(), "PI/2 is more than 1!");
    assert!(
        half_pi - half > T::from_u8(1).unwrap(),
        "PI/2 is more than 1.5!"
    );
    if int_bits >= 2 {
        assert!(half_pi < T::from_u8(2).unwrap(), "PI/2 is less than 2!");
    }
    if N < 49 {
        assert_eq!(utils::raw_u64(&half_pi), pi_raw_u64_48 >> (49 - N));
    } else {
        assert_eq!(utils::raw_u64(&half_pi) >> (N - 49), pi_raw_u64_48);
    }

    if int_bits >= 2 {
        let r_pi = T::FRAC_1_PI();
        let pi = T::PI();
        let r = r_pi * pi;
        assert!(utils::nearly_equal(r, T::ONE, 2), "1/PI * PI should be 1");

        let r = half_pi * T::FRAC_2_PI();
        assert!(utils::nearly_equal(r, T::ONE, 2), "2/PI * PI/2 should be 1");
    }

    if int_bits >= 3 {
        let r_pi = T::FRAC_1_PI();
        let r = r_pi * T::TAU() * half;
        assert!(
            utils::nearly_equal(r, T::ONE, 2),
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

    if N < 48 {
        assert_eq!(utils::raw_u64(&T::E()), e_raw_u64_48 >> (48 - N));
        assert_eq!(utils::raw_u64(&T::LN_2()), ln2_raw_u64_48 >> (48 - N));
        assert_eq!(utils::raw_u64(&T::LN_10()), ln10_raw_u64_48 >> (48 - N));
    } else {
        assert_eq!(utils::raw_u64(&T::E()) >> (N - 48), e_raw_u64_48);
        assert_eq!(utils::raw_u64(&T::LN_2()) >> (N - 48), ln2_raw_u64_48);
        assert_eq!(utils::raw_u64(&T::LN_10()) >> (N - 48), ln10_raw_u64_48 - 1);
        // ah hack the rounding for now
    }

    assert!(
        utils::nearly_equal(T::LOG2_E() * T::LN_2(), T::ONE, 2),
        "logs"
    );
    assert!(
        utils::nearly_equal(T::LOG2_10() * T::LN_2(), T::LN_10(), 2),
        "logs"
    );

    assert!(
        utils::nearly_equal(T::LOG10_E() * T::LN_10(), T::ONE, 2),
        "logs"
    );
    assert!(
        utils::nearly_equal(T::LOG10_2() * T::LN_10(), T::LN_2(), 2),
        "logs"
    );
}
