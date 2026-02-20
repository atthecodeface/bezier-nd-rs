use thiserror::Error;

/// An error from the Bezier library
#[derive(Error, Debug)]
pub enum BezierError {
    /// A matrix has been provided to a function that has an incorrect size
    #[error("bad matrix size provided - needed {0}, but was provided {1}")]
    BadMatrixSize(usize, usize),
    /// A vector has been provided to a function that has an incorrect length
    #[error("bad vector size provided - needed {0}, but was provided {1}")]
    BadVectorSize(usize, usize),
    /// Bezier build has been given a set of constraints that either are not
    /// linearly independent (and hence are underconstained) or are
    /// overconstraining (e.g. more constraints provided than control points in
    /// the eventual Bezier)
    #[error("bezier build was either underconstrained or overconstrained")]
    BadBuildConstraints,
    /// A Bezier was attempted to be built, but it would have required a Bezier
    /// of a higher degree than the type itself supports
    #[error("bezier build would make a degree {0} Bezier but the most possible was {1}")]
    MaxBuildDegree(usize, usize),
    /// A Bezier build must provide enough coefficients (from constraints) to
    /// create the control points for the Bezier given its degree, but there
    /// were too few
    #[error("bezier build needed {0} coefficients (as it has one fewer degree) but was provided with too few coefficients {1}")]
    BuildTooFewCoefficients(usize, usize),
    /// A Build constraint was provided that is not supported currently by the
    /// library, or was just invalid (a fifth derivative cannot be requested for
    /// a linear Bezier, for exanple)
    #[error("bezier build was provided with a constraint that could not be built (derivative of more than 1 degree is not implemented yet)")]
    UnableToBuild(),
}
