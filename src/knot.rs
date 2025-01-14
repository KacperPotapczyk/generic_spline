use std::{collections::HashMap, error::Error, fmt::Display};

pub struct Knot {
    pub x: f64,
    pub y: f64,
    pub continuity: usize,
    pub derivatives_values: HashMap<usize, f64>
}

impl Knot {
    pub fn new(x: f64, y: f64, continuity: usize, derivatives_values: HashMap<usize, f64>) -> Result<Self, Box<dyn Error>> {

        let max_derivative_degree = derivatives_values.iter()
                .map(|kv| *kv.0)
                .max()
                .unwrap_or(0);
        if continuity < max_derivative_degree {
            return Err(Box::new(
                KnotError("maximal derivatives_values degree is greater than continuity".to_string())
            ))
        }

        Ok(Knot { x, y, continuity, derivatives_values })
    }

    pub fn c0(x: f64, y: f64) -> Self {
        Knot { x, y, continuity: 0, derivatives_values: HashMap::new() }
    }

    pub fn c1(x: f64, y: f64) -> Self {
        Knot { x, y, continuity: 1, derivatives_values: HashMap::new() }
    }

    pub fn c2(x: f64, y: f64) -> Self {
        Knot { x, y, continuity: 2, derivatives_values: HashMap::new() }
    }

    pub fn cn(x: f64, y:f64, continuity: usize) -> Self {
        Knot { x, y, continuity, derivatives_values: HashMap::new() }
    }

    pub fn fix0(x: f64, y: f64) -> Self {
        Knot { x, y, continuity: 0, derivatives_values: HashMap::new() }
    }

    pub fn fix1(x: f64, y: f64, first_derivative: f64) -> Self {
        Knot { x, y, continuity: 1, derivatives_values: HashMap::from([(1, first_derivative)]) }
    }

    pub fn fix2(x: f64, y: f64, first_derivative: f64, second_derivative: f64) -> Self {
        Knot { x, y, continuity: 2, derivatives_values: HashMap::from([(1, first_derivative), (2, second_derivative)]) }
    }

    pub fn fixn(x: f64, y: f64, derivatives_values: HashMap<usize, f64>) -> Self {
        let max_derivative_degree = derivatives_values.iter()
                .map(|kv| *kv.0)
                .max()
                .unwrap_or(0);
        Knot { x, y, continuity: max_derivative_degree, derivatives_values }
    }

}

impl Ord for Knot {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.x.total_cmp(&other.x)
    }
}

impl PartialOrd for Knot {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Knot {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x
    }
}

impl Eq for Knot { }

#[derive(Debug)]
struct KnotError(String);

impl Display for KnotError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Error in Knot: {}", self.0)
    }
}

impl Error for KnotError {}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let x = 1.0;
        let y = 2.5;
        let continuity = 2;
        let derivatives_values = HashMap::from([(1, 4.0)]);
        let knot = Knot::new(x, y, continuity, derivatives_values).unwrap();

        assert_eq!(x, knot.x);
        assert_eq!(y, knot.y);
        assert_eq!(continuity, knot.continuity);
        assert_eq!(1, knot.derivatives_values.len());
        assert_eq!(4.0, *knot.derivatives_values.get(&1).unwrap());
    }

    #[test]
    fn test_c0() {
        let x = 1.0;
        let y = 2.5;
        let knot = Knot::c0(x, y);

        assert_eq!(x, knot.x);
        assert_eq!(y, knot.y);
        assert_eq!(0, knot.continuity);
        assert_eq!(0, knot.derivatives_values.len());
    }

    #[test]
    fn test_c1() {
        let x = 1.0;
        let y = 2.5;
        let knot = Knot::c1(x, y);

        assert_eq!(x, knot.x);
        assert_eq!(y, knot.y);
        assert_eq!(1, knot.continuity);
        assert_eq!(0, knot.derivatives_values.len());
    }

    #[test]
    fn test_c2() {
        let x = 1.0;
        let y = 2.5;
        let knot = Knot::c2(x, y);

        assert_eq!(x, knot.x);
        assert_eq!(y, knot.y);
        assert_eq!(2, knot.continuity);
        assert_eq!(0, knot.derivatives_values.len());
    }

    #[test]
    fn test_cn() {
        let x = 1.0;
        let y = 2.5;
        let continuity = 3;
        let knot = Knot::cn(x, y, continuity);

        assert_eq!(x, knot.x);
        assert_eq!(y, knot.y);
        assert_eq!(continuity, knot.continuity);
        assert_eq!(0, knot.derivatives_values.len());
    }

    #[test]
    fn test_fix0() {
        let x = 1.0;
        let y = 2.5;
        let knot = Knot::fix0(x, y);

        assert_eq!(x, knot.x);
        assert_eq!(y, knot.y);
        assert_eq!(0, knot.continuity);
        assert_eq!(0, knot.derivatives_values.len());
    }

    #[test]
    fn test_fix1() {
        let x = 1.0;
        let y = 2.5;
        let first_derivative = 0.5;
        let knot = Knot::fix1(x, y, first_derivative);

        assert_eq!(x, knot.x);
        assert_eq!(y, knot.y);
        assert_eq!(1, knot.continuity);
        assert_eq!(1, knot.derivatives_values.len());
        assert_eq!(first_derivative, *knot.derivatives_values.get(&1).unwrap());
    }

    #[test]
    fn test_fix2() {
        let x = 1.0;
        let y = 2.5;
        let first_derivative = 0.5;
        let second_derivative = -3.5;
        let knot = Knot::fix2(x, y, first_derivative, second_derivative);

        assert_eq!(x, knot.x);
        assert_eq!(y, knot.y);
        assert_eq!(2, knot.continuity);
        assert_eq!(2, knot.derivatives_values.len());
        assert_eq!(first_derivative, *knot.derivatives_values.get(&1).unwrap());
        assert_eq!(second_derivative, *knot.derivatives_values.get(&2).unwrap());
    }

    #[test]
    fn test_fixn() {
        let x = 1.0;
        let y = 2.5;
        let first_derivative = 0.5;
        let second_derivative = -3.5;
        let third_derivative = 5.0;

        let derivatives_values = HashMap::from([(1, first_derivative), (2, second_derivative), (3, third_derivative)]);

        let knot = Knot::fixn(x, y, derivatives_values);

        assert_eq!(x, knot.x);
        assert_eq!(y, knot.y);
        assert_eq!(3, knot.continuity);
        assert_eq!(3, knot.derivatives_values.len());
        assert_eq!(first_derivative, *knot.derivatives_values.get(&1).unwrap());
        assert_eq!(second_derivative, *knot.derivatives_values.get(&2).unwrap());
        assert_eq!(third_derivative, *knot.derivatives_values.get(&3).unwrap());
    }
}