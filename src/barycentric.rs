use crate::*;

use faer::traits::ComplexField;
use num_complex::{ComplexFloat};

impl<T> Barycentric<T>
where T: ComplexField + ComplexFloat
{
    pub fn new(nodes: Vec<T>, values: Vec<T>, weights: Vec<T>, error: f64) -> Self {
        assert_eq!(nodes.len(), values.len());
        assert_eq!(nodes.len(), weights.len());
        let weightsxvals = weights.iter().zip(values.iter()).map(|(w,v)| *w * *v).collect();
        Barycentric { nodes, values, weights, weightsxvals, error }
    }
}