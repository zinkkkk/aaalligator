use faer::traits::ComplexField;
use num_complex::ComplexFloat;
use num_traits::NumAssignOps;
use crate::{Barycentric, continuum::aaa_continuum};

pub trait BaryFromFn<T>
{
    fn into_bary(self) -> Barycentric<T>;
}

impl<T, F> BaryFromFn<T> for F
where T: ComplexFloat + ComplexField + NumAssignOps<T>,
F: Fn(T) -> T + std::marker::Sync + Copy, f64: From<<T as ComplexField>::Real>,
{
    fn into_bary(self) -> Barycentric<T> {
    aaa_continuum(self, None, 130, 0, 1e-13, 130, 2)
    }
}
