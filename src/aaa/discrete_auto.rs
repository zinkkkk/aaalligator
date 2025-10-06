use core::f64;
use iter_num_tools::lin_space;
use crate::{Barycentric, discrete::aaa_discreet};

pub fn auto_discrete_f_range(f: Vec<f64>) -> Barycentric<f64> {

    let min = f.iter().cloned().reduce(f64::min).unwrap();
    let max = f.iter().cloned().reduce(f64::max).unwrap();

    let flen = f.len();
    let pts: Vec<f64> = lin_space(min..=max, flen).collect();
    aaa_discreet(f, pts, 50, 1000.0*f64::EPSILON, 10)
}

pub fn auto_discrete_f_zero_one(f: Vec<f64>) -> Barycentric<f64> {
    let flen = f.len() as i32;
    let pts: Vec<f64> = lin_space(0.0..=1.0, flen as usize).collect();
    aaa_discreet(f, pts, ((flen / 2) - 1).try_into().unwrap(), f64::EPSILON, 10)
}

pub fn auto_discrete_f_negone_one(f: Vec<f64>) -> Barycentric<f64> {
    let flen = f.len() as i32;
    let pts: Vec<f64> = lin_space(-1.0..=1.0, flen as usize).collect();
    aaa_discreet(f, pts, ((flen / 2) - 1).try_into().unwrap(), 1e-8, (flen / 2) - 1)
}