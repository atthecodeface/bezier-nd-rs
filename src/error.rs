use thiserror::Error;

#[derive(Error, Debug)]
pub enum BezierError {
    #[error("bad matrix size provided - needed {0}, but was provided {1}")]
    BadMatrixSize(usize, usize),
    #[error("bad vector size provided - needed {0}, but was provided {1}")]
    BadVectorSize(usize, usize),
    #[error("bezier build was either underconstrained or overconstrained")]
    BadBuildConstraints,
    #[error("bezier build would make a degree {0} Bezier but the most possible was {1}")]
    MaxBuildDegree(usize, usize),
    #[error("bezier build needed {0} coefficients (as it has one fewer degree) but was provided with too few coefficients {1}")]
    BuildTooFewCoefficients(usize, usize),
    #[error("bezier build was provided with a constraint that could not be built (derivative of more than 1 degree is not implemented yet)")]
    UnableToBuild(),
    #[error("unknown bezier error")]
    Unknown,
}
