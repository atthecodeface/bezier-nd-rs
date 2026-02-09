mod coeffs;
mod elevate_matrix;

/// Functions to permit creation of cubic Bezier curves of circular arcs
pub mod arc;

/// Functions that split Bernestein Bezier curves, generating new control points
pub mod de_casteljau;

/// Functions to calculate point values and first derivatives, for &[[F;D]]
pub mod values;

/// Functions for distances
pub mod distance;

pub use coeffs::*;
pub use elevate_matrix::*;
