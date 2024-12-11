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
                .max();
        if continuity < max_derivative_degree.unwrap_or(0) {
            return Err(Box::new(
                KnotError("maximal derivatives_values degree is greater than continuity".to_string())
            ))
        }

        Ok(Knot { x, y, continuity, derivatives_values })
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