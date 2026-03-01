#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ArithCode {
    /// Indicates a perfectly valid result was obtained
    Ok,
    /// Indicates result exceeded the range of the most positive value; the actual result provided is ideally the Wrapping value
    ///
    /// If saturating, then set to the most positive value; if checked, return None; if wrapping, then ignore; if plain, then panic
    ///
    /// For values that support inifinites, then the result may be replaced with +Infinity
    OverflowMax,
    /// Indicates result exceeded the range of the most negative value; the actual result provided is ideally the Wrapping value
    ///
    /// If saturating, then set to the most negative value; if checked, return None; if wrapping, then ignore; if plain, then panic
    ///
    /// For values that support inifinites, then the result may be replaced with +Infinity
    OverflowMin,
    /// A division by zero
    ///
    /// If saturating, then set to the most negative value; if checked, return None; if wrapping, then ignore; if plain, then panic
    DivideByZero,
    /// Not a number
    ///
    /// This may be returned for example for a square root of a negative number, or asin of 2, or log of a negative number, etc
    NotANumber,
}
