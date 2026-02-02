mod coeffs;
mod elevate_matrix;

/// Functions to evaluate various metrics for Bernstein polynomial Bezier curves
pub mod metrics;

/// Functions that split Bernestein Bezier curves, generating new control points
pub mod split;
pub use coeffs::*;
pub use elevate_matrix::*;
