extern crate generic_spline;

use generic_spline::{Knot, Spline};

fn main() {

    let x_min = 0.0;
    let x_max = 6.0;

    let knots = vec![
        Knot::fix0(x_min, 1.0),
        Knot::c0(1.0, -1.0),
        Knot::c0(2.0, 0.0),
        Knot::c0(4.0, 3.0),
        Knot::c0(5.0, 1.0),
        Knot::fix0(x_max, 1.0)
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