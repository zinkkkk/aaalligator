use crate::*;

use iter_num_tools::lin_space;
use faer::traits::ComplexField;
use num_complex::ComplexFloat;
use num_traits::NumAssignOps;

pub trait BaryFromFn<T>
{
    /// simple conversion of functions or vecs into barycentric approimations
    fn into_bary(self) -> Barycentric<T>;
}

impl<T, F> BaryFromFn<T> for F
where T: ComplexFloat + ComplexField + NumAssignOps<T>,
F: Fn(T) -> T + std::marker::Sync + Copy, f64: From<<T as ComplexField>::Real>,
{
    fn into_bary(self) -> Barycentric<T> {
    aaa_continuum(self, None, 150, 0, 1e-13, 150, 3)
    }
}

impl<T> BaryFromFn<T> for Vec<T>
where T: ComplexFloat + ComplexField + NumAssignOps<T> + num_traits::FromPrimitive, f64: From<<T as ComplexField>::Real>
{
    fn into_bary(self) -> Barycentric<T> {
    let len = self.len();
    let z: Vec<T> = lin_space(T::from(-1.0).unwrap()..=T::from(1.0).unwrap(), len).collect();
    aaa_discreet(self, z, len, 1e-13, len)
    }
}

impl<T> BaryFromFn<T> for (Vec<T>, Vec<T>)
where T: ComplexFloat + ComplexField + NumAssignOps<T> + num_traits::FromPrimitive, f64: From<<T as ComplexField>::Real>
{
    fn into_bary(self) -> Barycentric<T> {
    assert!(self.0.len() == self.1.len());
    let len = self.0.len();
    aaa_discreet(self.0, self.1, len, 1e-13, len)
    }
}