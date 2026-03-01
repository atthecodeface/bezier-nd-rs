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
