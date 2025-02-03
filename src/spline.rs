use std::{collections::HashMap, error::Error, fmt::Display};

use nalgebra::{DMatrix, DVector};

use crate::{knot::Knot, polynomial::Polynomial};

pub struct Spline {
    knots: Vec<Knot>,
    polynomials: Vec<Polynomial>,
    min_x: f64,
    max_x: f64,
    is_knot_spacing_uniform: bool,
    last_interval_cache: usize,
}

impl Spline {
    pub fn new(knots: Vec<Knot>) -> Result<Self, Box<dyn Error>> {

        if knots.len() < 2 {
            return Err(Box::new(SplineError("Spline must have at least 2 knots".to_string())));
        }

        let number_of_intervals = knots.len() - 1;
        let mut spline = Spline {
            knots,
            polynomials: Vec::with_capacity(number_of_intervals),
            min_x: 0.0,
            max_x: 0.0,
            is_knot_spacing_uniform: false,
            last_interval_cache: 0,
        };

        spline.sort_knots();
        spline.check_knots_spacing()?;
        spline.calculate_polynomials()?;
        return Ok(spline);
    }

    pub fn interpolate(&self, x: f64) -> Result<f64, Box<dyn Error>> {
        if self.is_in_range(x) {
            let index = self.find_interval_index(x);
            Ok(self.polynomials[index].evaluate(x))
        } else {
            return Err(Box::new(SplineError("x is out of range".to_string())));
        }
    }

    pub fn cached_interpolate(&mut self, x: f64) -> Result<f64, Box<dyn Error>> {
        if self.is_in_range(x) {
            let index = self.find_interval_index_with_cache(x);
            Ok(self.polynomials[index].evaluate(x))
        } else {
            return Err(Box::new(SplineError("x is out of range".to_string())));
        }
    }

    pub fn batch_interpolate(&self, x_vector: &Vec<f64>) -> Result<Vec<f64>, Box<dyn Error>> {

        if x_vector.iter().all(|x| self.is_in_range(*x)) {

            let mut results = Vec::with_capacity(x_vector.len());
            let mut index = 0;

            for x in x_vector {
                index = self.find_interval_index_with_hint(index, *x);
                results.push(self.polynomials[index].evaluate(*x));
            }
            return Ok(results);
        } else {
            return Err(Box::new(SplineError("x is out of range".to_string())));
        }
    }

    pub fn extrapolate(&self, x: f64) -> f64 {

        match self.evaluate_on_boundaries(x) {
            Some(result) => return result,
            None => {
                let index = self.find_interval_index(x);
                return self.polynomials[index].evaluate(x)
            },
        }
    }

    pub fn cached_extrapolate(&mut self, x: f64) -> f64 {

        match self.evaluate_on_boundaries(x) {
            Some(result) => return result,
            None => {
                let index = self.find_interval_index_with_cache(x);
                return self.polynomials[index].evaluate(x)
            },
        }
    }

    pub fn batch_extrapolate(&self, x_vector: &Vec<f64>) -> Vec<f64> {

        let mut results = Vec::with_capacity(x_vector.len());
        let mut index = 0;

        for x in x_vector {
            match self.evaluate_on_boundaries(*x) {
                Some(result) => results.push(result),
                None => {
                    index = self.find_interval_index_with_hint(index, *x);
                    results.push(self.polynomials[index].evaluate(*x));
                },
            }
        }
        return results;
    }

    fn sort_knots(&mut self) {
        self.knots.sort();
        self.min_x = self.knots[0].x;
        self.max_x = self.knots[self.knots.len() - 1].x;
    }

    fn check_knots_spacing(&mut self) -> Result<(), Box<dyn Error>> {

        let x_spacing_vec: Vec<f64> = self.knots.iter()
            .map(|k| k.x)
            .collect::<Vec<f64>>()
            .windows(2)
            .map(|w| w[1] - w[0])
            .collect();

        if x_spacing_vec.iter().any(|spacing| *spacing < 1e-16) {
            return Err(Box::new(SplineError("Knots have equal x values".to_string())));
        }

        self.is_knot_spacing_uniform = x_spacing_vec
            .windows(2)
            .map(|spacing| (spacing[1] - spacing[0]).abs())
            .all(|difference| difference < 1e-16);
    
        Ok(())
    }
    
    fn calculate_polynomials(&mut self) -> Result<(), Box<dyn Error>> {
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

        let polynomial_derivative_coefficients =
            self.calculate_polynomial_derivative_coefficients(&intervals_degree);

        let mut x1_pow: Vec<f64> = (0..intervals_coefficients_number[0])
            .into_iter()
            .map(|p| self.knots[0].x.powi(p as i32))
            .collect();
        let mut x0_pow ;

        // println!("intervals_degree: {:?}", intervals_degree);

        for i in 0..number_of_intervals {
            let index_start = intervals_index_start[i];
            let number_of_coefficients = intervals_coefficients_number[i];

            if x1_pow.len() >= number_of_coefficients {
                x0_pow = x1_pow.clone();
            } else {
                x0_pow = (0..number_of_coefficients)
                .into_iter()
                .map(|p| self.knots[i].x.powi(p as i32))
                .collect();
            }
            
            x1_pow = (0..number_of_coefficients)
                .into_iter()
                .map(|p| self.knots[i + 1].x.powi(p as i32))
                .collect();

            self.function_value_equation_coefficients(
                number_of_coefficients,
                index_start,
                &x0_pow,
                self.knots[i].y,
                &mut equation_counter,
                &mut matrix,
                &mut rhs,
            );

            self.function_value_equation_coefficients(
                number_of_coefficients,
                index_start,
                &x1_pow,
                self.knots[i + 1].y,
                &mut equation_counter,
                &mut matrix,
                &mut rhs,
            );

            if i < number_of_intervals - 1 {
                self.continuity_equation_coefficients(
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
                    &polynomial_derivative_coefficients,
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
                &polynomial_derivative_coefficients,
            );

            self.derivatives_values_equation_coefficients(
                number_of_coefficients,
                index_start,
                &x1_pow,
                &self.knots[i + 1].derivatives_values,
                &mut equation_counter,
                &mut matrix,
                &mut rhs,
                &polynomial_derivative_coefficients,
            );
        }

        // println!("matrix: {}", matrix);
        // println!("rhs: {}", rhs);

        let solution = match matrix.lu().solve(&rhs) {
            Some(solution) => solution,
            None => return Err(Box::new(SplineError("Error while solving set of equations".to_string()))),
        };
        // println!("solution: {}", solution);

        for i in 0..number_of_intervals {
            let number_of_coefficients = intervals_coefficients_number[i];
            let interval_index_start = intervals_index_start[i];
            self.create_polynomial_for_interval(
                number_of_coefficients,
                interval_index_start,
                &solution,
            );
        }
        Ok(())
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

    fn calculate_polynomial_derivative_coefficients(
        &mut self,
        intervals_degree: &Vec<usize>,
    ) -> DMatrix<f64> {
        let max_interval_degree = *intervals_degree.iter().max().unwrap_or(&0);
        let max_derivative = self.knots.iter().map(|k| k.continuity).max().unwrap_or(0);

        let mut polynomial_derivative_coefficents =
            DMatrix::<f64>::zeros(max_interval_degree + 1, max_derivative + 1);
        for i in 0..=max_interval_degree {
            for d in 0..=max_derivative {
                polynomial_derivative_coefficents[(i, d)] =
                    self.polynomial_derivative_coefficient(i, d);
            }
        }

        polynomial_derivative_coefficents
    }

    fn polynomial_derivative_coefficient(
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

    fn function_value_equation_coefficients(
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

    fn continuity_equation_coefficients(
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
        polynomial_derivative_coefficients: &DMatrix<f64>,
    ) {
        for order in 1..=continiuity {
            if derivatives_values.contains_key(&order) {
                continue;
            }

            for c in 0..number_of_coefficients_0 {
                matrix[(*equation_counter, index_start_0 + c)] = self
                    .derivative_equation_coefficient(
                        c,
                        order,
                        x_pow,
                        polynomial_derivative_coefficients,
                    );
            }

            for c in 0..number_of_coefficients_1 {
                matrix[(*equation_counter, index_start_1 + c)] = -1.0
                    * self.derivative_equation_coefficient(
                        c,
                        order,
                        x_pow,
                        polynomial_derivative_coefficients,
                    );
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
        polynomial_derivative_coefficients: &DMatrix<f64>,
    ) {
        for (order, value) in derivatives_values {
            for c in 0..number_of_coefficients {
                matrix[(*equation_counter, index_start + c)] = self
                    .derivative_equation_coefficient(
                        c,
                        *order,
                        &x_pow,
                        polynomial_derivative_coefficients,
                    );
            }
            rhs[*equation_counter] = *value;
            *equation_counter += 1;
        }
    }

    fn derivative_equation_coefficient(
        &self,
        polynomial_order: usize,
        derievative_order: usize,
        x_pow: &Vec<f64>,
        polynomial_derivative_coefficients: &DMatrix<f64>,
    ) -> f64 {
        if polynomial_order < derievative_order {
            return 0.0;
        } else {
            return x_pow[polynomial_order - derievative_order]
                * polynomial_derivative_coefficients[(polynomial_order, derievative_order)];
        }
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
        if self.is_knot_spacing_uniform {
            return self.find_interval_index_uniform(x);
        } else {
            return self.find_interval_index_bisect(x);
        }
    }

    fn find_interval_index_bisect(&self, x:f64) -> usize {
        let size = self.knots.len();
        let mut min = 0;
        let mut max = size - 1;
        
        while max - min > 1 {
            let mid = (min + max) / 2;
            if x < self.knots[mid].x {
                max = mid;
            } else {
                min = mid;
            }
        }
        return min;
    }

    fn find_interval_index_uniform(&self, x:f64) -> usize {

        let relative_x = (x - self.min_x) / (self.max_x - self.min_x);
        if relative_x < 0.0 || relative_x > 1.0 {
            panic!("x is out of range")
        }
        let index = (relative_x * (self.knots.len() - 1) as f64).floor() as usize;
        if index == self.knots.len() - 1 {
            return index - 1;
        }
        return index;
    }

    fn find_interval_index_with_cache(&mut self, x: f64) -> usize {
        
        if !self.is_in_interval_range(self.last_interval_cache, x) {

            if self.last_interval_cache < self.knots.len()-1 && self.is_in_interval_range(self.last_interval_cache+1, x) {
                self.last_interval_cache += 1;
            } else {
                self.last_interval_cache = self.find_interval_index(x);
            }
        }
        self.last_interval_cache
    }

    fn find_interval_index_with_hint(&self, index_hint: usize, x: f64) -> usize {
        
        if !self.is_in_interval_range(index_hint, x) {

            if index_hint < self.knots.len()-1 && self.is_in_interval_range(index_hint+1, x) {
                return index_hint + 1;
            } else {
                return self.find_interval_index(x);
            }
        }
        return index_hint;
    }

    fn is_in_interval_range(&self, interval_index: usize, x: f64) -> bool {
        self.knots[interval_index].x <= x && x <= self.knots[interval_index+1].x
    }

    fn evaluate_on_boundaries(&self, x:f64) -> Option<f64> {
        let size = self.knots.len();
        if x < self.knots[1].x {
            Some(self.polynomials[0].evaluate(x))
        } else if x > self.knots[size-2].x {
            Some(self.polynomials[size-2].evaluate(x))
        } else {
            None
        }
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

        let knots = vec![
            Knot::new(0.0, 0.0_f64.powi(2), 1, HashMap::from([(1, 0.0)])).unwrap(),
            Knot::new(0.9, 0.9_f64.powi(2), 2, HashMap::new()).unwrap(),
            Knot::new(1.1, 1.1_f64.powi(2), 2, HashMap::new()).unwrap(),
            Knot::new(1.7, 1.7_f64.powi(2), 2, HashMap::new()).unwrap(),
            Knot::new(2.0, 2.0_f64.powi(2), 1, HashMap::from([(1, 4.0)])).unwrap()
        ];

        let spline = Spline::new(knots).unwrap();

        assert!(!spline.is_knot_spacing_uniform);

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

        assert!(spline.interpolate(3.0).is_err());

        assert_approx_eq!(spline.extrapolate(-1.0), (-1.0_f64).powi(2), eps);
        assert_approx_eq!(spline.extrapolate(-0.2), (-0.2_f64).powi(2), eps);
        assert_approx_eq!(spline.extrapolate(2.7), 2.7_f64.powi(2), eps);
        assert_approx_eq!(spline.extrapolate(3.0), 3.0_f64.powi(2), eps);
        assert_approx_eq!(spline.extrapolate(1.5), 1.5_f64.powi(2), eps);
    }

    #[test]
    fn over_x_squared_function_with_cache() {
        // knots lay on f(x) = x^2 function
        let eps = 1e-6;

        let knots = vec![
            Knot::new(0.0, 0.0_f64.powi(2), 1, HashMap::from([(1, 0.0)])).unwrap(),
            Knot::new(0.9, 0.9_f64.powi(2), 2, HashMap::new()).unwrap(),
            Knot::new(1.1, 1.1_f64.powi(2), 2, HashMap::new()).unwrap(),
            Knot::new(1.7, 1.7_f64.powi(2), 2, HashMap::new()).unwrap(),
            Knot::new(2.0, 2.0_f64.powi(2), 1, HashMap::from([(1, 4.0)])).unwrap()
        ];

        let mut spline = Spline::new(knots).unwrap();

        assert!(!spline.is_knot_spacing_uniform);

        assert_approx_eq!(spline.cached_interpolate(0.0).unwrap(), 0.0, eps);
        assert_approx_eq!(spline.cached_interpolate(0.13).unwrap(), 0.13_f64.powi(2), eps);
        assert_approx_eq!(spline.cached_interpolate(0.69).unwrap(), 0.69_f64.powi(2), eps);
        assert_approx_eq!(spline.cached_interpolate(1.0).unwrap(), 1.0, eps);
        assert_approx_eq!(spline.cached_interpolate(1.13).unwrap(), 1.13_f64.powi(2), eps);
        assert_approx_eq!(
            spline.cached_interpolate(1.8643128).unwrap(),
            1.8643128_f64.powi(2),
            eps
        );
        assert_approx_eq!(spline.cached_interpolate(2.0).unwrap(), 4.0, eps);

        assert!(spline.cached_interpolate(3.0).is_err());

        assert_approx_eq!(spline.cached_extrapolate(0.0001), (0.0001_f64).powi(2), eps);
        assert_approx_eq!(spline.cached_extrapolate(0.0002), (0.0002_f64).powi(2), eps);
        assert_approx_eq!(spline.cached_extrapolate(1.11), (1.11_f64).powi(2), eps);
        assert_approx_eq!(spline.cached_extrapolate(-1.0), (-1.0_f64).powi(2), eps);
        assert_approx_eq!(spline.cached_extrapolate(-0.2), (-0.2_f64).powi(2), eps);
        assert_approx_eq!(spline.cached_extrapolate(2.7), 2.7_f64.powi(2), eps);
        assert_approx_eq!(spline.cached_extrapolate(3.0), 3.0_f64.powi(2), eps);
        assert_approx_eq!(spline.cached_extrapolate(1.5), 1.5_f64.powi(2), eps);
    }



    #[test]
    fn over_x_squared_function_batch() {
        // knots lay on f(x) = x^2 function
        let eps = 1e-6;

        let knots = vec![
            Knot::new(0.0, 0.0_f64.powi(2), 1, HashMap::from([(1, 0.0)])).unwrap(),
            Knot::new(0.9, 0.9_f64.powi(2), 2, HashMap::new()).unwrap(),
            Knot::new(1.1, 1.1_f64.powi(2), 2, HashMap::new()).unwrap(),
            Knot::new(1.7, 1.7_f64.powi(2), 2, HashMap::new()).unwrap(),
            Knot::new(2.0, 2.0_f64.powi(2), 1, HashMap::from([(1, 4.0)])).unwrap()
        ];

        let spline = Spline::new(knots).unwrap();

        assert!(!spline.is_knot_spacing_uniform);

        let x_vector = vec![0.0, 0.13, 0.69, 1.0, 1.13, 1.8643128, 2.0];
        let result = spline.batch_interpolate(&x_vector).unwrap();

        assert_eq!(x_vector.len(), result.len());
        for i in 0..x_vector.len() {
            assert_approx_eq!(result[i], x_vector[i].powi(2), eps);
        }

        let x_vector = vec![0.0, 0.13, 0.69, 1.0, 3.0];
        assert!(spline.batch_interpolate(&x_vector).is_err());

        let x_vector = vec![0.0001, 0.0002, 1.11, -1.0, -0.2, 2.7, 3.0];
        let result = spline.batch_extrapolate(&x_vector);

        assert_eq!(x_vector.len(), result.len());
        for i in 0..x_vector.len() {
            assert_approx_eq!(result[i], x_vector[i].powi(2), eps);
        }
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

        let spline = Spline::new(knots).unwrap();

        assert!(spline.is_knot_spacing_uniform);

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

        let spline = Spline::new(knots).unwrap();

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

        let spline = Spline::new(knots).unwrap();

        assert_approx_eq!(spline.interpolate(0.0).unwrap(), y0, eps);
        assert_approx_eq!(spline.interpolate(0.2).unwrap(), 4.12, eps);
        assert_approx_eq!(spline.interpolate(0.5).unwrap(), 3.25, eps);
        assert_approx_eq!(spline.interpolate(1.0).unwrap(), y1, eps);
        assert_approx_eq!(spline.interpolate(1.5).unwrap(), 2.75, eps);
        assert_approx_eq!(spline.interpolate(1.8).unwrap(), 1.88, eps);
        assert_approx_eq!(spline.interpolate(2.0).unwrap(), y2, eps);

        assert_approx_eq!(spline.extrapolate(-1.0), 13.0, eps);
        assert_approx_eq!(spline.extrapolate(-0.2), 6.12, eps);
        assert_approx_eq!(spline.extrapolate(2.7), -3.97, eps);
        assert_approx_eq!(spline.extrapolate(3.0), -7.0, eps);
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

        let spline = Spline::new(knots).unwrap();

        assert_approx_eq!(spline.interpolate(0.0).unwrap(), y0, eps);
        assert_approx_eq!(spline.interpolate(0.2).unwrap(), 4.08, eps);
        assert_approx_eq!(spline.interpolate(0.5).unwrap(), 3.0, eps);
        assert_approx_eq!(spline.interpolate(1.0).unwrap(), y1, eps);
        assert_approx_eq!(spline.interpolate(1.5).unwrap(), 2.0, eps);
        assert_approx_eq!(spline.interpolate(1.8).unwrap(), 1.52, eps);
        assert_approx_eq!(spline.interpolate(2.0).unwrap(), y2, eps);
    }

    #[test]
    fn tree_point_c0_c0_c2_fix() {
        let eps = 1e-6;
        let y0 = -1.0;
        let y1 = 3.0;
        let y2 = -1.0;

        let knot1 = Knot::new(0.0, y0, 0, HashMap::new()).unwrap();
        let knot2 = Knot::new(1.0, y1, 0, HashMap::new()).unwrap();
        let knot3 = Knot::new(2.0, y2, 2, HashMap::from([(1, -10.0), (2, -16.0)])).unwrap();
        let knots = vec![knot1, knot2, knot3];

        let spline = Spline::new(knots).unwrap();

        assert_approx_eq!(spline.interpolate(0.0).unwrap(), y0, eps);
        assert_approx_eq!(spline.interpolate(0.2).unwrap(), -0.2, eps);
        assert_approx_eq!(spline.interpolate(0.5).unwrap(), 1.0, eps);
        assert_approx_eq!(spline.interpolate(1.0).unwrap(), y1, eps);
        assert_approx_eq!(spline.interpolate(1.5).unwrap(), 2.25, eps);
        assert_approx_eq!(spline.interpolate(1.8).unwrap(), 0.696, eps);
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

        let spline = Spline::new(knots).unwrap();

        assert_approx_eq!(spline.interpolate(0.0).unwrap(), y0, eps);
        assert_approx_eq!(spline.interpolate(0.2).unwrap(), 4.48, eps);
        assert_approx_eq!(spline.interpolate(0.5).unwrap(), 4.0, eps);
        assert_approx_eq!(spline.interpolate(1.0).unwrap(), y1, eps);
        assert_approx_eq!(spline.interpolate(1.5).unwrap(), 5.0, eps);
        assert_approx_eq!(spline.interpolate(1.8).unwrap(), 6.08, eps);
        assert_approx_eq!(spline.interpolate(2.0).unwrap(), y2, eps);
    }

    #[test]
    fn two_point_c1_fix_c1_fix() {
        let eps = 1e-6;
        let y0 = -4.0;
        let y1 = 6.0;

        let knot1 = Knot::new(0.0, y0, 1, HashMap::from([(1, 3.0)])).unwrap();
        let knot2 = Knot::new(2.0, y1, 1, HashMap::from([(1, 1.0)])).unwrap();
        let knots = vec![knot1, knot2];

        let spline = Spline::new(knots).unwrap();

        assert!(spline.is_knot_spacing_uniform);

        assert_approx_eq!(spline.interpolate(0.0).unwrap(), y0, eps);
        assert_approx_eq!(spline.interpolate(0.2).unwrap(), -3.252, eps);
        assert_approx_eq!(spline.interpolate(0.5).unwrap(), -1.6875, eps);
        assert_approx_eq!(spline.interpolate(1.5).unwrap(), 4.4375, eps);
        assert_approx_eq!(spline.interpolate(2.0).unwrap(), y1, eps);


        assert_approx_eq!(spline.extrapolate(-1.0), -1.5, eps);
        assert_approx_eq!(spline.extrapolate(-0.2), -4.428, eps);
        assert_approx_eq!(spline.extrapolate(2.7), 3.7355, eps);
        assert_approx_eq!(spline.extrapolate(3.0), 0.5, eps);
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

        let spline = Spline::new(knots).unwrap();

        assert!(spline.is_knot_spacing_uniform);

        assert_approx_eq!(spline.interpolate(0.0).unwrap(), y0, eps);
        assert_approx_eq!(spline.interpolate(0.2).unwrap(), 2.344, eps);
        assert_approx_eq!(spline.interpolate(0.5).unwrap(), 1.375, eps);
        assert_approx_eq!(spline.interpolate(1.0).unwrap(), y1, eps);
        assert_approx_eq!(spline.interpolate(1.5).unwrap(), 3.0, eps);
        assert_approx_eq!(spline.interpolate(1.8).unwrap(), 4.008, eps);
        assert_approx_eq!(spline.interpolate(2.0).unwrap(), y2, eps);

        assert!(spline.interpolate(3.0).is_err());

        assert_approx_eq!(spline.extrapolate(-1.0), 1.0, eps);
        assert_approx_eq!(spline.extrapolate(-0.2), 3.496, eps);
        assert_approx_eq!(spline.extrapolate(2.7), -4.848, eps);
        assert_approx_eq!(spline.extrapolate(3.0), -15.0, eps);
    }

    #[test]
    fn tree_point_c1_fix_c2_c1_fix_with_cache() {
        let eps = 1e-6;
        let y0 = 3.0;
        let y1 = 1.0;
        let y2 = 4.0;

        let knot1 = Knot::new(0.0, y0, 1, HashMap::from([(1, -3.0)])).unwrap();
        let knot2 = Knot::new(1.0, y1, 2, HashMap::new()).unwrap();
        let knot3 = Knot::new(2.0, y2, 1, HashMap::from([(1, -2.0)])).unwrap();
        let knots = vec![knot1, knot2, knot3];

        let mut spline = Spline::new(knots).unwrap();

        assert!(spline.is_knot_spacing_uniform);

        assert_approx_eq!(spline.cached_interpolate(0.0).unwrap(), y0, eps);
        assert_approx_eq!(spline.cached_interpolate(0.2).unwrap(), 2.344, eps);
        assert_approx_eq!(spline.cached_interpolate(0.5).unwrap(), 1.375, eps);
        assert_approx_eq!(spline.cached_interpolate(1.0).unwrap(), y1, eps);
        assert_approx_eq!(spline.cached_interpolate(1.5).unwrap(), 3.0, eps);
        assert_approx_eq!(spline.cached_interpolate(1.8).unwrap(), 4.008, eps);
        assert_approx_eq!(spline.cached_interpolate(2.0).unwrap(), y2, eps);
        
        assert!(spline.cached_interpolate(3.0).is_err());

        assert_approx_eq!(spline.cached_extrapolate(-1.0), 1.0, eps);
        assert_approx_eq!(spline.cached_extrapolate(-0.2), 3.496, eps);
        assert_approx_eq!(spline.cached_extrapolate(2.7), -4.848, eps);
        assert_approx_eq!(spline.cached_extrapolate(3.0), -15.0, eps);
    }

    #[test]
    fn tree_point_c1_fix_c2_c1_fix_batch() {
        let eps = 1e-6;
        let y0 = 3.0;
        let y1 = 1.0;
        let y2 = 4.0;

        let knot1 = Knot::new(0.0, y0, 1, HashMap::from([(1, -3.0)])).unwrap();
        let knot2 = Knot::new(1.0, y1, 2, HashMap::new()).unwrap();
        let knot3 = Knot::new(2.0, y2, 1, HashMap::from([(1, -2.0)])).unwrap();
        let knots = vec![knot1, knot2, knot3];

        let spline = Spline::new(knots).unwrap();

        assert!(spline.is_knot_spacing_uniform);

        let x_vector = vec![0.0, 0.2, 0.5, 1.0, 1.5, 1.8, 2.0];
        let expected_result = vec![y0, 2.344, 1.375, y1, 3.0, 4.008, y2];
        let result = spline.batch_interpolate(&x_vector).unwrap();

        assert_eq!(x_vector.len(), result.len());
        for i in 0..x_vector.len() {
            assert_approx_eq!(result[i], expected_result[i], eps);
        }

        let x_vector = vec![-1.0, -0.2, 2.7, 3.0];
        let expected_result = vec![1.0, 3.496, -4.848, -15.0];
        let result = spline.batch_extrapolate(&x_vector);

        assert_eq!(x_vector.len(), result.len());
        for i in 0..x_vector.len() {
            assert_approx_eq!(result[i], expected_result[i], eps);
        }
    }

    #[test]
    fn tree_point_c1_fix_c2_c1_fix_non_uniform() {
        let eps = 1e-6;
        let y0 = -4.0;
        let y1 = 6.0;
        let y2 = -4.0;

        let knot1 = Knot::new(0.0, y0, 1, HashMap::from([(1, 3.0)])).unwrap();
        let knot2 = Knot::new(2.0, y1, 2, HashMap::new()).unwrap();
        let knot3 = Knot::new(3.0, y2, 1, HashMap::from([(1, -27.0)])).unwrap();
        let knots = vec![knot1, knot2, knot3];

        let spline = Spline::new(knots).unwrap();

        assert!(!spline.is_knot_spacing_uniform);

        assert_approx_eq!(spline.interpolate(0.0).unwrap(), y0, eps);
        assert_approx_eq!(spline.interpolate(0.2).unwrap(), -3.252, eps);
        assert_approx_eq!(spline.interpolate(0.5).unwrap(), -1.6875, eps);
        assert_approx_eq!(spline.interpolate(1.5).unwrap(), 4.4375, eps);
        assert_approx_eq!(spline.interpolate(2.0).unwrap(), y1, eps);
        assert_approx_eq!(spline.interpolate(2.3).unwrap(), 5.688, eps);
        assert_approx_eq!(spline.interpolate(2.8).unwrap(), 0.528, eps);
        assert_approx_eq!(spline.interpolate(3.0).unwrap(), y2, eps);

        assert_approx_eq!(spline.extrapolate(-1.0), -1.5, eps);
        assert_approx_eq!(spline.extrapolate(-0.2), -4.428, eps);
        assert_approx_eq!(spline.extrapolate(3.2), -10.368, eps);
        assert_approx_eq!(spline.extrapolate(4.0), -60.0, eps);
        assert_approx_eq!(spline.extrapolate(0.2), -3.252, eps);
        assert_approx_eq!(spline.extrapolate(2.3), 5.688, eps);
        assert_approx_eq!(spline.extrapolate(2.8), 0.528, eps);
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

        let spline = Spline::new(knots).unwrap();

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

        let spline = Spline::new(knots).unwrap();

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

        let spline = Spline::new(knots).unwrap();

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

        let spline = Spline::new(knots).unwrap();

        assert_approx_eq!(spline.interpolate(0.0).unwrap(), y0, eps);
        assert_approx_eq!(spline.interpolate(0.2).unwrap(), 2.4272, eps);
        assert_approx_eq!(spline.interpolate(0.5).unwrap(), 2.375, eps);
        assert_approx_eq!(spline.interpolate(1.0).unwrap(), y1, eps);
        assert_approx_eq!(spline.interpolate(1.5).unwrap(), 0.625, eps);
        assert_approx_eq!(spline.interpolate(1.8).unwrap(), 2.1328, eps);
        assert_approx_eq!(spline.interpolate(2.0).unwrap(), y2, eps);

        assert_approx_eq!(spline.extrapolate(-1.0), -1.0, eps);
        assert_approx_eq!(spline.extrapolate(-0.2), 1.2592, eps);
        assert_approx_eq!(spline.extrapolate(2.7), 14.4538, eps);
        assert_approx_eq!(spline.extrapolate(3.0), 19.0, eps);
    }

    #[test]
    fn test_one_knot_error() {
        let knot1 = Knot::c0(0.0, 2.0);
        let knots = vec![knot1];

        let spline = Spline::new(knots);

        assert!(spline.is_err())
    }

    #[test]
    fn test_equal_x_knot_values() {
        let y0 = 2.0;
        let y1 = 1.0;
        let y2 = 4.0;

        let knot1 = Knot::c0(0.0, y0);
        let knot2 = Knot::c0(0.0, y1);
        let knot3 = Knot::c0(1.0, y2);
        let knots = vec![knot1, knot2, knot3];

        let spline = Spline::new(knots);

        assert!(spline.is_err())
    }

    #[test]
    fn example() {

        let x_min = 0.0;
        let x_max = 6.0;

        let knots = vec![
            Knot::fix1(x_min, 1.0, 0.0),
            Knot::c2(1.0, -1.0),
            Knot::c2(2.0, 0.0),
            Knot::c2(3.0, -1.0),
            Knot::c2(4.0, 3.0),
            Knot::c2(5.0, 0.5),
            Knot::fix1(x_max, 1.0, 0.0)
        ];

        let spline = Spline::new(knots).unwrap();

        let number_of_points = 60;
        let step = (x_max - x_min) /  number_of_points as f64;
        for i in 0..=number_of_points {
            let x = x_min + step*i as f64;
            println!("{:.2};{:.2}", x, spline.interpolate(x).unwrap());
        }
        assert!(true);
    }

    #[ignore]
    #[test]
    fn perfomance() {
        use std::time::Instant;
        use rand::Rng;

        let x_min = 0.0;
        let x_max = 6.0;
        let mut rng = rand::thread_rng();

        let mut knots = vec![
            Knot::fix1(x_min, 1.0, 0.0),
            Knot::fix1(x_max, 1.0, 0.0)
        ];

        let knots_number = 30;
        let knot_step = (x_max - x_min) / knots_number as f64;

        for i in 1..knots_number {
            let x = x_min + knot_step*i as f64;
            let y = rng.gen_range(0.0..10.0);
            knots.push(Knot::c2(x, y));
        }

        let mut spline = Spline::new(knots).unwrap();

        let number_of_points = 300;
        let step = (x_max - x_min) / number_of_points as f64;

        let mut x_vector = Vec::new();
        for i in 10..=number_of_points {
            x_vector.push(x_min + step*i as f64);
        }

        let now = Instant::now();
        for x in x_vector.iter() {
            assert!(spline.interpolate(*x).unwrap() >= -10.0);
        }
        let elapsed = now.elapsed();
        println!("interpolate time: {:.2?}", elapsed);

        let now = Instant::now();
        for x in x_vector.iter() {
            assert!(spline.cached_interpolate(*x).unwrap() >= -10.0);
        }
        let elapsed = now.elapsed();
        println!("cached_interpolate time: {:.2?}", elapsed);

        let now = Instant::now();
        let result = spline.batch_interpolate(&x_vector).unwrap();
        assert!(result.len() == x_vector.len());
        let elapsed = now.elapsed();
        println!("batch_interpolate time: {:.2?}", elapsed);
    }
}
