extern crate generic_spline;

use generic_spline::{Knot, Spline};

fn main() {

    let knots = vec![
        Knot::fix1(1.0, 1.0, 0.0),
        Knot::c2(2.0, 0.0),
        Knot::c2(3.0, -2.0),
        Knot::c2(4.0, 1.0),
        Knot::fix1(5.0, 1.0, -1.0)
    ];

    let spline = Spline::new(knots).unwrap();

    let x_min = 0.0;
    let x_max = 6.0;
    let number_of_steps = 60;
    let step = (x_max - x_min) / number_of_steps as f64;

    println!("x;y");
    for i in 0..=number_of_steps {
        let x = x_min + step * i as f64;
        println!("{:.2};{:.2}", x, spline.extrapolate(x));
    }    
}