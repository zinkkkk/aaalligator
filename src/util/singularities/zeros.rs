use faer::{diag::Diag, matrix_free::LinOp, prelude::*, traits::{RealField}};
use num_complex::{Complex, ComplexFloat};
use num_traits::{float::FloatCore};

use crate::Barycentric;

pub trait Zeros<T> {
    fn zeros(&self) -> Vec<Complex<T>>;
}

impl<T> Zeros<T> for Barycentric<T>
where T: RealField + ComplexFloat + FloatCore
{
    /// finds zeros from an approximation
    fn zeros(&self) -> Vec<Complex<T>> {
        let m: usize = self.weights.len();
        let mut e: Mat<T> = Mat::zeros(m + 1, m + 1);

        for i in 0..m {
            e[(i + 1, 0)] = T::from(1.0).unwrap();
            e[(0, i + 1)] = self.weightsxvals[i];
            e[(i + 1, i + 1)] = self.nodes[i];
        }

        let mut b: Mat<T> = Mat::identity(m + 1, m + 1);
        b[(0, 0)] = T::from(0.0).unwrap();

        let gevd = e.generalized_eigen(b).unwrap();

        let s_a = gevd.S_a();
        let s_b = gevd.S_b();

        let mut eigs = Diag::zeros(s_a.nrows());
        zip!(&mut eigs, &s_a, &s_b).for_each(|unzip!(a, b, c)| *a = b / *c);

        let eigs = eigs.column_vector();

        let zer = eigs
            .iter()
            .cloned()
            .filter(|&f | f.is_finite())
            .collect::<Vec<Complex<T>>>();

        zer
    }
}