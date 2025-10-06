use faer::traits::ComplexField;
use num_complex::Complex;
use num_traits::{float::FloatCore, zero};

use crate::*;

pub trait Residues<T> {
    fn residues(&self, poles: &[Complex<T>]) -> Vec<Complex<T>>;
}

impl<T> Residues<T> for Barycentric<T>
where T: ComplexField + FloatCore
{
    fn residues(&self, poles: &[Complex<T>]) -> Vec<Complex<T>> {

        let zj = &self.nodes;
        let fj = &self.values;
        let wj = &self.weights;

        let n = |t|
            zj
                .iter()
                .zip(fj.iter())
                .zip(wj.iter())
                .map(|((x, y), w)| Complex::from(*y * *w) / (t - Complex::new(*x, zero())))
                .sum::<Complex<T>>();

        let ddiff = |t|
            zj
                .iter()
                .zip(wj.iter())
                .map(|(&x, &w)| Complex::from(w.neg()) / t - Complex::new(x, zero()).powi(2))
                .sum::<Complex<T>>();

        let res: Vec<Complex<T>> = poles
            .iter()
            .map(|&p| n(p) / ddiff(p))
            .collect();

        res
    }
}