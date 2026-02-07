mod coeffs;
mod elevate_matrix;

/// Functions to permit creation of cubic Bezier curves of circular arcs
pub mod arc;

/// Functions that split Bernestein Bezier curves, generating new control points
pub mod de_casteljau;
pub use coeffs::*;
pub use elevate_matrix::*;
