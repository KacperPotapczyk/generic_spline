use std::{collections::HashMap, error::Error, fmt::Display};

use nalgebra::{DMatrix, DVector};

use crate::{knot::Knot, polynomial::Polynomial};

pub struct Spline {
    knots: Vec<Knot>,
    polynomials: Vec<Polynomial>,
    min_x: f64,
    max_x: f64,
}

impl Spline {
    pub fn new(knots: Vec<Knot>) -> Self {
        let number_of_intervals = knots.len() - 1;
        let mut spline = Spline {
            knots,
            polynomials: Vec::with_capacity(number_of_intervals),
            min_x: 0.0,
            max_x: 0.0,
        };

        spline.sort_knots();
        spline.calculate_polynomials();
        return spline;
    }

    pub fn interpolate(&self, x: f64) -> Result<f64, Box<dyn Error>> {
        if self.is_in_range(x) {
            let index = self.find_interval_index(x);
            let result = self.polynomials[index].evaluate(x);
            Ok(result)
        } else {
            return Err(Box::new(SplineError("x is out of range".to_string())));
        }
    }

    fn sort_knots(&mut self) {
        self.knots.sort();
        self.min_x = self.knots[0].x;
        self.max_x = self.knots[self.knots.len() - 1].x;
    }

    fn calculate_polynomials(&mut self) {
        let number_of_intervals = self.knots.len() - 1;

        let intervals_degree = self.calculate_intervals_degree(number_of_intervals);

        let intervals_coefficients_number: Vec<usize> =
            intervals_degree.iter().map(|d| *d + 1).collect();

        let matrix_size = intervals_coefficients_number.iter().sum();
        let mut matrix = DMatrix::<f64>::zeros(matrix_size, matrix_size);
        let mut rhs = DVector::<f64>::zeros(matrix_size);

        let mut equation_counter = 0;
        let intervals_index_start = self
            .calculate_intervals_index_start(number_of_intervals, &intervals_coefficients_number);

        let mut x1_pow: Vec<f64> = (0..intervals_coefficients_number[0])
            .into_iter()
            .map(|p| self.knots[0].x.powi(p as i32))
            .collect();

        for i in 0..number_of_intervals {
            let index_start = intervals_index_start[i];
            let number_of_coefficients = intervals_coefficients_number[i];

            let x0_pow: Vec<f64> = x1_pow.clone();
            x1_pow = (0..number_of_coefficients)
                .into_iter()
                .map(|p| self.knots[i + 1].x.powi(p as i32))
                .collect();

            self.function_value_equation_coefficents(
                number_of_coefficients,
                index_start,
                &x0_pow,
                self.knots[i].y,
                &mut equation_counter,
                &mut matrix,
                &mut rhs,
            );

            self.function_value_equation_coefficents(
                number_of_coefficients,
                index_start,
                &x1_pow,
                self.knots[i + 1].y,
                &mut equation_counter,
                &mut matrix,
                &mut rhs,
            );

            if i < number_of_intervals - 1 {
                self.continuity_equation_coefficents(
                    self.knots[i + 1].continuity,
                    number_of_coefficients,
                    intervals_coefficients_number[i + 1],
                    index_start,
                    intervals_index_start[i + 1],
                    &x1_pow,
                    &self.knots[i + 1].derivatives_values,
                    &mut equation_counter,
                    &mut matrix,
                    &mut rhs,
                );
            }

            self.derivatives_values_equation_coefficients(
                number_of_coefficients,
                index_start,
                &x0_pow,
                &self.knots[i].derivatives_values,
                &mut equation_counter,
                &mut matrix,
                &mut rhs,
            );

            self.derivatives_values_equation_coefficients(
                number_of_coefficients,
                index_start,
                &x1_pow,
                &self.knots[i + 1].derivatives_values,
                &mut equation_counter,
                &mut matrix,
                &mut rhs,
            );
        }

        println!("intervals_degree: {:?}", intervals_degree);
        println!("matrix: {}", matrix);
        println!("rhs: {}", rhs);

        let solution = matrix.lu().solve(&rhs).unwrap();
        println!("solution: {}", solution);

        for i in 0..number_of_intervals {
            let number_of_coefficients = intervals_coefficients_number[i];
            let interval_index_start = intervals_index_start[i];
            self.create_polynomial_for_interval(
                number_of_coefficients,
                interval_index_start,
                &solution,
            );
        }
    }

    fn calculate_intervals_degree(&mut self, number_of_intervals: usize) -> Vec<usize> {
        let mut intervals_degree: Vec<usize> = Vec::with_capacity(number_of_intervals);
        for i in 0..number_of_intervals {
            let mut degree: usize = 1;

            degree += self.knots[i].derivatives_values.len();
            degree += self.knots[i + 1].derivatives_values.len();

            intervals_degree.push(degree);
        }

        // correct degree using continuity without fixed derivative
        for i in 0..number_of_intervals - 1 {
            let continuity_minus_fixed_values =
                self.knots[i + 1].continuity - self.knots[i + 1].derivatives_values.len();
            let degree_from_continuity = continuity_minus_fixed_values / 2;

            if continuity_minus_fixed_values % 2 == 0 {
                intervals_degree[i] += degree_from_continuity;
                intervals_degree[i + 1] += degree_from_continuity;
            } else {
                if intervals_degree[i] > intervals_degree[i + 1] {
                    intervals_degree[i] += degree_from_continuity;
                    intervals_degree[i + 1] += degree_from_continuity + 1;
                } else {
                    intervals_degree[i] += degree_from_continuity + 1;
                    intervals_degree[i + 1] += degree_from_continuity;
                }
            }
        }
        intervals_degree
    }

    fn calculate_intervals_index_start(
        &self,
        number_of_intervals: usize,
        intervals_coefficients_number: &Vec<usize>,
    ) -> Vec<usize> {
        let mut intervals_index_start: Vec<usize> = Vec::with_capacity(number_of_intervals);

        let mut index_start = 0;
        for interval_coefficient_number in intervals_coefficients_number.iter() {
            intervals_index_start.push(index_start);
            index_start += interval_coefficient_number;
        }
        intervals_index_start
    }

    fn function_value_equation_coefficents(
        &self,
        number_of_coefficients: usize,
        index_start: usize,
        x_pow: &Vec<f64>,
        y_value: f64,
        equation_counter: &mut usize,
        matrix: &mut DMatrix<f64>,
        rhs: &mut DVector<f64>,
    ) {
        for c in 0..number_of_coefficients {
            matrix[(*equation_counter, index_start + c)] = x_pow[c];
        }
        rhs[*equation_counter] = y_value;
        *equation_counter += 1;
    }

    fn continuity_equation_coefficents(
        &self,
        continiuity: usize,
        number_of_coefficients_0: usize,
        number_of_coefficients_1: usize,
        index_start_0: usize,
        index_start_1: usize,
        x_pow: &Vec<f64>,
        derivatives_values: &HashMap<usize, f64>,
        equation_counter: &mut usize,
        matrix: &mut DMatrix<f64>,
        rhs: &mut DVector<f64>,
    ) {
        for order in 1..=continiuity {
            if derivatives_values.contains_key(&order) {
                continue;
            }

            for c in 0..number_of_coefficients_0 {
                matrix[(*equation_counter, index_start_0 + c)] =
                    self.derivative_equation_coefficent(c, order, x_pow);
            }

            for c in 0..number_of_coefficients_1 {
                matrix[(*equation_counter, index_start_1 + c)] =
                    -1.0 * self.derivative_equation_coefficent(c, order, x_pow);
            }

            rhs[*equation_counter] = 0.0;
            *equation_counter += 1;
        }
    }

    fn derivatives_values_equation_coefficients(
        &self,
        number_of_coefficients: usize,
        index_start: usize,
        x_pow: &Vec<f64>,
        derivatives_values: &HashMap<usize, f64>,
        equation_counter: &mut usize,
        matrix: &mut DMatrix<f64>,
        rhs: &mut DVector<f64>,
    ) {
        for (order, value) in derivatives_values {
            for c in 0..number_of_coefficients {
                matrix[(*equation_counter, index_start + c)] =
                    self.derivative_equation_coefficent(c, *order, &x_pow);
            }
            rhs[*equation_counter] = *value;
            *equation_counter += 1;
        }
    }

    fn derivative_equation_coefficent(
        &self,
        polynomial_order: usize,
        derievative_order: usize,
        x_pow: &Vec<f64>,
    ) -> f64 {
        if polynomial_order < derievative_order {
            return 0.0;
        } else {
            return x_pow[polynomial_order - derievative_order]
                * self.polynomial_derivative_coeff(polynomial_order, derievative_order);
        }
    }

    fn polynomial_derivative_coeff(
        &self,
        polynomial_order: usize,
        derievative_order: usize,
    ) -> f64 {
        let mut multiplier = 1.0;
        let mut coeff = polynomial_order as f64;
        for _ in 0..derievative_order {
            multiplier *= coeff;
            coeff -= 1.0;
        }
        multiplier
    }

    fn create_polynomial_for_interval(
        &mut self,
        number_of_coefficients: usize,
        interval_index_start: usize,
        solution: &DVector<f64>,
    ) {
        let mut coefficients = Vec::with_capacity(number_of_coefficients);

        for c in interval_index_start..interval_index_start + number_of_coefficients {
            coefficients.push(solution[c]);
        }

        let interval_polynomial = Polynomial::new(coefficients);
        self.polynomials.push(interval_polynomial);
    }

    fn is_in_range(&self, x: f64) -> bool {
        self.min_x <= x && x <= self.max_x
    }

    fn find_interval_index(&self, x: f64) -> usize {
        let size = self.knots.len();
        for i in 1..size - 1 {
            if x < self.knots[i].x {
                return i - 1;
            }
        }
        return size - 2;
    }
}

#[derive(Debug)]
struct SplineError(String);

impl Display for SplineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Error in Spline: {}", self.0)
    }
}

impl Error for SplineError {}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;
    use std::collections::HashMap;

    use super::*;

    #[test]
    fn over_x_squared_function() {
        // knots lay on f(x) = x^2 function
        let eps = 1e-6;

        let knot1 = Knot::new(0.0, 0.0_f64.powi(2), 1, HashMap::from([(1, 0.0)])).unwrap();
        let knot2 = Knot::new(1.0, 1.0_f64.powi(2), 2, HashMap::new()).unwrap();
        let knot3 = Knot::new(2.0, 2.0_f64.powi(2), 1, HashMap::from([(1, 4.0)])).unwrap();
        let knots = vec![knot1, knot2, knot3];

        let spline = Spline::new(knots);

        assert_approx_eq!(spline.interpolate(0.0).unwrap(), 0.0, eps);
        assert_approx_eq!(spline.interpolate(0.13).unwrap(), 0.13_f64.powi(2), eps);
        assert_approx_eq!(spline.interpolate(0.69).unwrap(), 0.69_f64.powi(2), eps);
        assert_approx_eq!(spline.interpolate(1.0).unwrap(), 1.0, eps);
        assert_approx_eq!(spline.interpolate(1.13).unwrap(), 1.13_f64.powi(2), eps);
        assert_approx_eq!(
            spline.interpolate(1.8643128).unwrap(),
            1.8643128_f64.powi(2),
            eps
        );
        assert_approx_eq!(spline.interpolate(2.0).unwrap(), 4.0, eps);
    }

    #[test]
    fn tree_point_c0_c0_c0() {
        let eps = 1e-6;
        let y0 = 4.0;
        let y1 = 2.0;
        let y2 = 6.0;

        let knot1 = Knot::new(0.0, y0, 0, HashMap::new()).unwrap();
        let knot2 = Knot::new(1.0, y1, 0, HashMap::new()).unwrap();
        let knot3 = Knot::new(2.0, y2, 0, HashMap::new()).unwrap();
        let knots = vec![knot1, knot2, knot3];

        let spline = Spline::new(knots);

        assert_approx_eq!(spline.interpolate(0.0).unwrap(), y0, eps);
        assert_approx_eq!(spline.interpolate(0.25).unwrap(), 3.5, eps);
        assert_approx_eq!(spline.interpolate(0.5).unwrap(), 3.0, eps);
        assert_approx_eq!(spline.interpolate(1.0).unwrap(), y1, eps);
        assert_approx_eq!(spline.interpolate(1.5).unwrap(), 4.0, eps);
        assert_approx_eq!(spline.interpolate(1.75).unwrap(), 5.0, eps);
        assert_approx_eq!(spline.interpolate(2.0).unwrap(), y2, eps);
    }

    #[test]
    fn tree_point_c0_c1_c0() {
        let eps = 1e-6;
        let y0 = 5.0;
        let y1 = 3.0;
        let y2 = 4.0;

        let knot1 = Knot::new(0.0, y0, 0, HashMap::new()).unwrap();
        let knot2 = Knot::new(1.0, y1, 1, HashMap::new()).unwrap();
        let knot3 = Knot::new(2.0, y2, 0, HashMap::new()).unwrap();
        let knots = vec![knot1, knot2, knot3];

        let spline = Spline::new(knots);

        assert_approx_eq!(spline.interpolate(0.0).unwrap(), y0, eps);
        assert_approx_eq!(spline.interpolate(0.2).unwrap(), 4.12, eps);
        assert_approx_eq!(spline.interpolate(0.5).unwrap(), 3.25, eps);
        assert_approx_eq!(spline.interpolate(1.0).unwrap(), y1, eps);
        assert_approx_eq!(spline.interpolate(1.5).unwrap(), 3.5, eps);
        assert_approx_eq!(spline.interpolate(1.8).unwrap(), 3.8, eps);
        assert_approx_eq!(spline.interpolate(2.0).unwrap(), y2, eps);
    }

    #[test]
    fn tree_point_c0_c1_fix_c0() {
        let eps = 1e-6;
        let y0 = 5.0;
        let y1 = 3.0;
        let y2 = 1.0;

        let knot1 = Knot::new(0.0, y0, 0, HashMap::new()).unwrap();
        let knot2 = Knot::new(1.0, y1, 1, HashMap::from([(1, 1.0)])).unwrap();
        let knot3 = Knot::new(2.0, y2, 0, HashMap::new()).unwrap();
        let knots = vec![knot1, knot2, knot3];

        let spline = Spline::new(knots);

        assert_approx_eq!(spline.interpolate(0.0).unwrap(), y0, eps);
        assert_approx_eq!(spline.interpolate(0.2).unwrap(), 4.12, eps);
        assert_approx_eq!(spline.interpolate(0.5).unwrap(), 3.25, eps);
        assert_approx_eq!(spline.interpolate(1.0).unwrap(), y1, eps);
        assert_approx_eq!(spline.interpolate(1.5).unwrap(), 2.75, eps);
        assert_approx_eq!(spline.interpolate(1.8).unwrap(), 1.88, eps);
        assert_approx_eq!(spline.interpolate(2.0).unwrap(), y2, eps);
    }

    #[test]
    fn tree_point_c1_fix_c0_c1_fix() {
        let eps = 1e-6;
        let y0 = 5.0;
        let y1 = 2.0;
        let y2 = 1.0;

        let knot1 = Knot::new(0.0, y0, 1, HashMap::from([(1, -5.0)])).unwrap();
        let knot2 = Knot::new(1.0, y1, 0, HashMap::new()).unwrap();
        let knot3 = Knot::new(2.0, y2, 1, HashMap::from([(1, -3.0)])).unwrap();
        let knots = vec![knot1, knot2, knot3];

        let spline = Spline::new(knots);

        assert_approx_eq!(spline.interpolate(0.0).unwrap(), y0, eps);
        assert_approx_eq!(spline.interpolate(0.2).unwrap(), 4.08, eps);
        assert_approx_eq!(spline.interpolate(0.5).unwrap(), 3.0, eps);
        assert_approx_eq!(spline.interpolate(1.0).unwrap(), y1, eps);
        assert_approx_eq!(spline.interpolate(1.5).unwrap(), 2.0, eps);
        assert_approx_eq!(spline.interpolate(1.8).unwrap(), 1.52, eps);
        assert_approx_eq!(spline.interpolate(2.0).unwrap(), y2, eps);
    }

    #[test]
    fn tree_point_c0_c2_c0() {
        let eps = 1e-6;
        let y0 = 5.0;
        let y1 = 4.0;
        let y2 = 7.0;

        let knot1 = Knot::new(0.0, y0, 0, HashMap::new()).unwrap();
        let knot2 = Knot::new(1.0, y1, 2, HashMap::new()).unwrap();
        let knot3 = Knot::new(2.0, y2, 0, HashMap::new()).unwrap();
        let knots = vec![knot1, knot2, knot3];

        let spline = Spline::new(knots);

        assert_approx_eq!(spline.interpolate(0.0).unwrap(), y0, eps);
        assert_approx_eq!(spline.interpolate(0.2).unwrap(), 4.48, eps);
        assert_approx_eq!(spline.interpolate(0.5).unwrap(), 4.0, eps);
        assert_approx_eq!(spline.interpolate(1.0).unwrap(), y1, eps);
        assert_approx_eq!(spline.interpolate(1.5).unwrap(), 5.0, eps);
        assert_approx_eq!(spline.interpolate(1.8).unwrap(), 6.08, eps);
        assert_approx_eq!(spline.interpolate(2.0).unwrap(), y2, eps);
    }

    #[test]
    fn tree_point_c1_fix_c2_c1_fix() {
        let eps = 1e-6;
        let y0 = 3.0;
        let y1 = 1.0;
        let y2 = 4.0;

        let knot1 = Knot::new(0.0, y0, 1, HashMap::from([(1, -3.0)])).unwrap();
        let knot2 = Knot::new(1.0, y1, 2, HashMap::new()).unwrap();
        let knot3 = Knot::new(2.0, y2, 1, HashMap::from([(1, -2.0)])).unwrap();
        let knots = vec![knot1, knot2, knot3];

        let spline = Spline::new(knots);

        assert_approx_eq!(spline.interpolate(0.0).unwrap(), y0, eps);
        assert_approx_eq!(spline.interpolate(0.2).unwrap(), 2.344, eps);
        assert_approx_eq!(spline.interpolate(0.5).unwrap(), 1.375, eps);
        assert_approx_eq!(spline.interpolate(1.0).unwrap(), y1, eps);
        assert_approx_eq!(spline.interpolate(1.5).unwrap(), 3.0, eps);
        assert_approx_eq!(spline.interpolate(1.8).unwrap(), 4.008, eps);
        assert_approx_eq!(spline.interpolate(2.0).unwrap(), y2, eps);
    }

    #[test]
    fn tree_point_c1_fix_c1_c1_fix() {
        let eps = 1e-6;
        let y0 = -3.0;
        let y1 = 6.0;
        let y2 = 13.0;

        let knot1 = Knot::new(0.0, y0, 1, HashMap::from([(1, 10.0)])).unwrap();
        let knot2 = Knot::new(1.0, y1, 1, HashMap::new()).unwrap();
        let knot3 = Knot::new(2.0, y2, 1, HashMap::from([(1, 9.0)])).unwrap();
        let knots = vec![knot1, knot2, knot3];

        let spline = Spline::new(knots);

        assert_approx_eq!(spline.interpolate(0.0).unwrap(), y0, eps);
        assert_approx_eq!(spline.interpolate(0.2).unwrap(), -0.944, eps);
        assert_approx_eq!(spline.interpolate(0.5).unwrap(), 2.125, eps);
        assert_approx_eq!(spline.interpolate(1.0).unwrap(), y1, eps);
        assert_approx_eq!(spline.interpolate(1.5).unwrap(), 9.0, eps);
        assert_approx_eq!(spline.interpolate(1.8).unwrap(), 11.28, eps);
        assert_approx_eq!(spline.interpolate(2.0).unwrap(), y2, eps);
    }

    #[test]
    fn tree_point_c1_fix_c1_fix_c1_fix() {
        let eps = 1e-6;
        let y0 = 4.0;
        let y1 = 3.0;
        let y2 = -7.0;

        let knot1 = Knot::new(0.0, y0, 1, HashMap::from([(1, 2.0)])).unwrap();
        let knot2 = Knot::new(1.0, y1, 1, HashMap::from([(1, -2.0)])).unwrap();
        let knot3 = Knot::new(2.0, y2, 1, HashMap::from([(1, -21.0)])).unwrap();
        let knots = vec![knot1, knot2, knot3];

        let spline = Spline::new(knots);

        assert_approx_eq!(spline.interpolate(0.0).unwrap(), y0, eps);
        assert_approx_eq!(spline.interpolate(0.2).unwrap(), 4.216, eps);
        assert_approx_eq!(spline.interpolate(0.5).unwrap(), 4.0, eps);
        assert_approx_eq!(spline.interpolate(1.0).unwrap(), y1, eps);
        assert_approx_eq!(spline.interpolate(1.5).unwrap(), 0.375, eps);
        assert_approx_eq!(spline.interpolate(1.8).unwrap(), -3.336, eps);
        assert_approx_eq!(spline.interpolate(2.0).unwrap(), y2, eps);
    }

    #[test]
    fn tree_point_c0_c2_fix2_c1_fix() {
        let eps = 1e-6;
        let y0 = -1.0;
        let y1 = 3.0;
        let y2 = -1.0;

        let knot1 = Knot::new(0.0, y0, 0, HashMap::new()).unwrap();
        let knot2 = Knot::new(1.0, y1, 2, HashMap::from([(2, -4.0)])).unwrap();
        let knot3 = Knot::new(2.0, y2, 1, HashMap::from([(1, -10.0)])).unwrap();
        let knots = vec![knot1, knot2, knot3];

        let spline = Spline::new(knots);

        assert_approx_eq!(spline.interpolate(0.0).unwrap(), y0, eps);
        assert_approx_eq!(spline.interpolate(0.2).unwrap(), 0.696, eps);
        assert_approx_eq!(spline.interpolate(0.5).unwrap(), 2.25, eps);
        assert_approx_eq!(spline.interpolate(1.0).unwrap(), y1, eps);
        assert_approx_eq!(spline.interpolate(1.5).unwrap(), 2.25, eps);
        assert_approx_eq!(spline.interpolate(1.8).unwrap(), 0.696, eps);
        assert_approx_eq!(spline.interpolate(2.0).unwrap(), y2, eps);
    }

    #[test]
    fn tree_point_c2_fix_c3_c1_fix() {
        let eps = 1e-6;
        let y0 = 2.0;
        let y1 = 1.0;
        let y2 = 4.0;

        let knot1 = Knot::new(0.0, y0, 2, HashMap::from([(1, 3.0), (2, -8.0)])).unwrap();
        let knot2 = Knot::new(1.0, y1, 3, HashMap::new()).unwrap();
        let knot3 = Knot::new(2.0, y2, 1, HashMap::from([(1, 11.0)])).unwrap();
        let knots = vec![knot1, knot2, knot3];

        let spline = Spline::new(knots);

        assert_approx_eq!(spline.interpolate(0.0).unwrap(), y0, eps);
        assert_approx_eq!(spline.interpolate(0.2).unwrap(), 2.4272, eps);
        assert_approx_eq!(spline.interpolate(0.5).unwrap(), 2.375, eps);
        assert_approx_eq!(spline.interpolate(1.0).unwrap(), y1, eps);
        assert_approx_eq!(spline.interpolate(1.5).unwrap(), 0.625, eps);
        assert_approx_eq!(spline.interpolate(1.8).unwrap(), 2.1328, eps);
        assert_approx_eq!(spline.interpolate(2.0).unwrap(), y2, eps);
    }
}
