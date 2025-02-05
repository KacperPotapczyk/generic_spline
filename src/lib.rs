//! Library of generic implementation of spline interpolation and extrapolation.
//! It does not assume that spline is used in graphics context.
//! 
//! # Example
//! ```
//! use generic_spline::{Knot, Spline};
//! use assert_approx_eq::assert_approx_eq;
//! 
//! let knots = vec![
//!     Knot::fix1(0.0, 3.0, -3.0),
//!     Knot::c2(1.0, 1.0),
//!     Knot::fix1(2.0, 4.0, -2.0)
//! ];
//! let spline = Spline::new(knots).unwrap();
//! 
//! assert_approx_eq!(2.344, spline.interpolate(0.2).unwrap(), 1e-6);
//! assert_approx_eq!(3.0, spline.interpolate(1.5).unwrap(), 1e-6);
//! ```

mod knot;
mod spline;
mod polynomial;

pub use knot::Knot;
pub use spline::Spline;