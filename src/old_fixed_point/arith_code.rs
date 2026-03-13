/// The code returned by some arithmentic operations to provide for saturating, wrapping, and checked results
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

impl ArithCode {
    /// Check that the code is okay, else panic with an appropriate message
    #[track_caller]
    pub fn check(self, reason: &str) {
        match self {
            ArithCode::Ok => (),
            ArithCode::NotANumber => {
                panic!("Not a number in {reason}");
            }
            ArithCode::OverflowMax => {
                panic!("Overflow (maximum) in {reason}");
            }
            ArithCode::OverflowMin => {
                panic!("Overflow (minimum) in {reason}");
            }
            ArithCode::DivideByZero => {
                panic!("Divide by zero in {reason}");
            }
        }
    }

    /// Check that the code is okay or overflowing, else panic with an appropriate message
    #[track_caller]
    pub fn wrapping(self, reason: &str) {
        match self {
            ArithCode::Ok => (),
            ArithCode::OverflowMax => (),
            ArithCode::OverflowMin => (),
            ArithCode::NotANumber => {
                panic!("Not a number in {reason}");
            }
            ArithCode::DivideByZero => {
                panic!("Divide by zero in {reason}");
            }
        }
    }

    /// Check that the code is okay or overflowing, else panic with an appropriate message
    #[track_caller]
    pub fn overflowing(self, reason: &str) -> bool {
        match self {
            ArithCode::Ok => false,
            ArithCode::OverflowMax => true,
            ArithCode::OverflowMin => true,
            ArithCode::NotANumber => {
                panic!("Not a number in {reason}");
            }
            ArithCode::DivideByZero => {
                panic!("Divide by zero in {reason}");
            }
        }
    }

    /// Check that the code is okay or overflowing, else panic with an appropriate message
    #[track_caller]
    pub fn checked<T>(self, value: T) -> Option<T> {
        match self {
            ArithCode::Ok => Some(value),
            _ => None,
        }
    }
}
