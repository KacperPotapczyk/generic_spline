extern crate generic_spline;

use generic_spline::{Knot, Spline};

fn main() {

    let x_min = 0.0;
    let x_max = 6.0;

    let knots = vec![
        Knot::fix1(x_min, 1.0, 0.0),
        Knot::c2(1.0, -1.0),
        Knot::c2(2.0, 0.0),
        Knot::c2(4.0, 3.0),
        Knot::c2(5.0, 1.0),
        Knot::fix1(x_max, 1.0, -1.0)
    ];

    let spline = Spline::new(knots).unwrap();

    let number_of_steps = 60;
    let step = (x_max - x_min) / number_of_steps as f64;

    let mut x_vector = Vec::new();

    for i in 0..=number_of_steps {
        x_vector.push(x_min + step * i as f64);
    }  

    let result = spline.batch_interpolate(&x_vector).unwrap();

    println!("x;y");
    for i in 0..=number_of_steps {
        println!("{:.2};{:.2}", x_vector[i], result[i]);
    }  
}