# generic_spline
Generic implementation of spline interpolation and extrapolation in Rust. It does not assume that spline is used in graphics context.
Spline is constructed over knots that can have constraints on required continuity and/or derivatives values. This implies that each spline interval (segment) can be polynomial of different order.

# Usage
Here we create spline with constraints on first derivative values on edge knots and C2 continuity on internal knots. After creating spline funcion value can be interpolated or extrapolated for given x.
``` rust
use generic_spline::{Knot, Spline};

let knots = vec![
    Knot::fix1(0.0, 1.0, 0.0),
    Knot::c2(2.0, 0.0),
    Knot::c2(5.0, 1.0),
    Knot::fix1(6.0, 1.0, -1.0)
];

let spline = Spline::new(knots).unwrap();

let result1 = spline.interpolate(3.5).unwrap();
let result2 = spline.interpolate(4.0).unwrap();
let result3 = spline.extrapolate(8.0);
```