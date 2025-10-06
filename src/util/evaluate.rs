use num_complex::ComplexFloat;
use num_traits::zero;

use crate::Barycentric;

pub trait Evaluate<T> {
    fn evaluate(&self, point: T) -> T;
}

pub trait EvaluateHandle<T> {
    fn evaluate_handle(&self) -> impl Fn(T) -> T;
}

impl<T> Evaluate<T> for Barycentric<T>
where
    T: num_traits::NumAssignOps
    + ComplexFloat
{
    fn evaluate(&self, point: T) -> T {

        if let Some(k) = self.nodes.iter().position(|&n| point == n) {
            return *self.values.get(k).unwrap();
        }

        let mut numerator = zero::<T>();
        let mut denominator = zero::<T>();

        for (n, w) in self.nodes.iter().zip(self.weights.iter()).zip(self.weightsxvals.iter()) {
        let ck = num_complex::ComplexFloat::recip(point - *n.0);
        numerator += ck * *w;
        denominator += ck * *n.1;
        }

    numerator / denominator
    }
}


impl<T> EvaluateHandle<T> for Barycentric<T>
where
    T: num_traits::NumAssignOps
    + ComplexFloat + std::iter::Sum
{
    /// returns a closure function to evaluate an approximation of type f(z)
    fn evaluate_handle(&self) -> impl Fn(T) -> T {
        |point: T| {


        if let Some(k) = self.nodes.iter().position(|&n| point == n) {
            return *self.values.get(k).unwrap();
        }

        let mut numerator = zero::<T>();
        let mut denominator = zero::<T>();

        for (n, w) in self.nodes.iter().zip(self.weights.iter()).zip(self.weightsxvals.iter()) {
        let ck = num_complex::ComplexFloat::recip(point - *n.0);
        numerator += ck * *w;
        denominator += ck * *n.1;
        }

    numerator / denominator
    }
    }
}