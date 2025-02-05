extern crate generic_spline;

use std::collections::HashMap;

use generic_spline::{Knot, Spline};

fn main() {

    let x_min = 0.0;
    let x_max = 6.0;

    let knots = vec![
        Knot::fix1(x_min, 1.0, 0.0),
        Knot::c2(1.0, -1.0),
        Knot::new(2.0, 0.0, 2, HashMap::from([(1, 0.0)])).unwrap(),
        Knot::c0(4.0, 3.0),
        Knot::c0(5.0, 1.0),
        Knot::fix2(x_max, 1.0, -1.0, -5.0)
    ];

    let spline = Spline::new(knots).unwrap();

    let number_of_steps = 60;
    let step = (x_max - x_min) / number_of_steps as f64;

    println!("x;y");
    for i in 0..=number_of_steps {
        let x = x_min + step * i as f64;
        println!("{:.2};{:.2}", x, spline.interpolate(x).unwrap());
    }
}