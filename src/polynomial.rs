pub struct Polynomial {
    coefficients: Vec<f64>,
    size: usize
}

impl Polynomial {

    pub fn new(coefficients: Vec<f64>) -> Self {
        let size = coefficients.len();
        Polynomial{ coefficients, size }
    }

    pub fn evaluate(&self, x: f64) -> f64 {
        let mut result = 0.0;
        for i in 0..self.size {
            result += x.powi(i as i32) * self.coefficients[i]
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;
    use super::*;

    #[test]
    fn evaluate() {

        let eps = 1e-6;
        let coefficients = vec![1.0, 2.5, -0.25];
        let polynomial = Polynomial::new(coefficients);

        assert_approx_eq!(polynomial.evaluate(2.1), 5.1475, eps);
        assert_approx_eq!(polynomial.evaluate(-3.14), -9.3149, eps);
        assert_approx_eq!(polynomial.evaluate(0.0), 1.0, eps);
    }
}