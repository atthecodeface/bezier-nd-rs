use super::Int;

/// This provides the constants required for trigonometry and logarithmic operations, for integer types
///
/// They are shifted such that they contain the maximum precision for their
/// type; so the value for `e` will (for a signed type) have a zero type bit
/// (sign bit) then a one bit (indicating two) then a zero units bit, etc
///
/// The following cannot be provided by the trait itself (as `>>` is not const)
/// but can be copied verbatim to implementations
///
pub trait UsefulConsts: Int {
    /// The value 'e' (2.718281828459045...) shifted to provide maximum precision
    const E: Self;

    /// The value `ln(2) == 1/log2(e)` (0.6931471805599453...) shifted to provide maximum precision
    ///
    /// `log2(x) = ln(x) / ln_2 = ln(x)*log2(e)`
    const LN_2: Self;
    /// The value `ln(10) == 1/log10(e)` (2.302585092994046...) shifted to provide maximum precision
    const LN_10: Self;

    /// The value `log2(e) == 1/ln(2)` (1.4426950408889634...) shifted to provide maximum precision
    ///
    /// `e^x == 2^(x.log2(e))`
    const LOG2_E: Self;
    /// The value `log2(10) == 1/log10(2)`  (3.3219280948873626...) shifted to provide maximum precision
    const LOG2_10: Self;

    /// The value `log10(2) == 1/log2(10)` (0.30102999566398114...) shifted to provide maximum precision
    const LOG10_2: Self;
    /// The value `log10(e) == 1/ln(10)` (0.43429448190325176...) shifted to provide maximum precision
    const LOG10_E: Self;

    /// The value 'sqrt(2)' (1.4142135623730951...) shifted to provide maximum precision
    ///
    /// Note: `1/sqrt(2) == sqrt(2)/2`
    const SQRT_2: Self;

    /// The value 'PI' (3.141592653589793...) shifted to provide maximum precision
    const PI: Self;
    /// The value '2/PI' (0.6366197723675814...) shifted to provide maximum precision
    const FRAC_2_PI: Self;
    /// The value 'PI/3' (1.0471975511965976...) shifted to provide maximum precision
    const FRAC_PI_3: Self;
    /// The value '1/sqrt(PI)' (0.5641895835477563...) shifted to provide maximum precision
    const FRAC_1_SQRT_PI: Self;
}

impl UsefulConsts for u128 {
    const E: Self = 0;
    const LN_2: Self = 0;
    const LN_10: Self = 0;
    const LOG2_E: Self = 0;
    const LOG2_10: Self = 0;
    const LOG10_2: Self = 0;
    const LOG10_E: Self = 0;
    // sqrt(2)/2 = 0xb504f333f9de6484597d89b3754abe9f1d6
    const SQRT_2: Self = 0xb504_f333_f9de_6484_597d_89b3_754a_be9f;
    // PI<<2 = 0x_c90fdaa22168c234c4c6628b80dc1cd129024e088a67cc74020bbea63b139b224
    const PI: Self = 0x_c90f_daa2_2168_c234_c4c6_628b_80dc_1cd1;
    // 1/PI = 0x_a2f9836e4e441529fc2757d1f534ddc0db6295993c439041fe5163abdebbc561c5
    const FRAC_2_PI: Self = 0x_a2f9_836e_4e44_1529_fc27_57d1_f534_ddc1;
    // PI/3 = 0x_860a91c16b9b2c232dd99707ab3d688b70ac3405b19a884d56b27f197cb7bcc18
    const FRAC_PI_3: Self = 0x_860a_91c1_6b9b_2c23_2dd9_9707_ab3d_688b;
    const FRAC_1_SQRT_PI: Self = 0;
}

macro_rules! useful_consts {
    {$t:ty, $ns:expr, $nb:expr} => {

    impl UsefulConsts for $t {
     const E: Self = (u128::E >> ($ns+128-$nb)) as $t;
     const LN_2: Self = (u128::LN_2 >> ($ns+128-$nb)) as $t;
     const LN_10: Self = (u128::LN_10 >> ($ns+128-$nb)) as $t;
     const LOG2_E: Self = (u128::LOG2_E >> ($ns+128-$nb)) as $t;
     const LOG2_10: Self = (u128::LOG2_10 >> ($ns+128-$nb)) as $t;
     const LOG10_2: Self = (u128::LOG10_2 >> ($ns+128-$nb)) as $t;
     const LOG10_E: Self = (u128::LOG10_E >> ($ns+128-$nb)) as $t;
     const SQRT_2: Self = (u128::SQRT_2 >> ($ns+128-$nb)) as $t;
     const PI: Self = (u128::PI >> ($ns+128-$nb)) as $t;
     const FRAC_2_PI: Self = (u128::FRAC_2_PI >> ($ns+128-$nb)) as $t;
     const FRAC_PI_3: Self = (u128::FRAC_PI_3 >> ($ns+128-$nb)) as $t;
     const FRAC_1_SQRT_PI: Self = (u128::FRAC_1_SQRT_PI >> ($ns+128-$nb)) as $t;
    }
    }
}

useful_consts!(i8, 1, 8);
useful_consts!(i16, 1, 16);
useful_consts!(i32, 1, 32);
useful_consts!(i64, 1, 64);
useful_consts!(i128, 1, 128);

useful_consts!(u8, 0, 8);
useful_consts!(u16, 0, 16);
useful_consts!(u32, 0, 32);
useful_consts!(u64, 0, 64);
