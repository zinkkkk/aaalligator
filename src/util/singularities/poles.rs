use crate::Barycentric;
use faer::{diag::Diag, matrix_free::LinOp, prelude::*, traits::{RealField, math_utils::zero}};
use num_complex::{Complex, ComplexFloat};
use num_traits::{float::FloatCore};

pub trait Poles<T> {
    fn poles(&self) -> Vec<Complex<T>>;
}

impl<T> Poles<T> for Barycentric<T>
where T: RealField + ComplexFloat + FloatCore
{
    /// finds poles from an approximation
    fn poles(&self) -> Vec<Complex<T>> {

        let m: usize = self.weights.len();
        let mut e: Mat<T> = Mat::zeros(m + 1, m + 1);

        for i in 0..m {
            e[(i + 1, 0)] = T::from(1.0).unwrap();
            e[(0, i + 1)] = self.weights[i];
            e[(i + 1, i + 1)] = self.nodes[i];
        }

        let mut b: Mat<T> = Mat::identity(m + 1, m + 1);
        b[(0, 0)] = zero();

        let gevd = e.generalized_eigen(b).unwrap();

        let s_a = gevd.S_a();
        let s_b = gevd.S_b();

        let mut eigs = Diag::zeros(s_a.nrows());
        zip!(&mut eigs, &s_a, &s_b).for_each(|unzip!(a, b, c)| *a = *b / *c);

        let eigs = eigs.column_vector();
    
        let pol = eigs
            .iter()
            .cloned()
            .filter(|&f | f.is_finite())
            .collect::<Vec<Complex<T>>>();

    pol
    }
}
