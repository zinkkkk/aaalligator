use faer::prelude::*;
use faer::traits::ComplexField;
use num_complex::{ComplexFloat};
use num_traits::{NumAssignOps, ToPrimitive, zero};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use crate::*;

/// AAA approximation from single in single out functions
pub fn aaa_continuum<T, F>(
    f: F,
    prenodes: Option<Vec<T>>,
    degree: i32,
    _lawson: i32,
    tol: f64,
    lookahead: i32,
    refinement: i32,
) -> Barycentric<T>
where T: ComplexFloat + ComplexField + NumAssignOps<T>,
F: Fn(T) -> T + std::marker::Sync, f64: From<<T as ComplexField>::Real>,
{

    let mut s: Vec<T> = vec![T::from(-1.0).unwrap(), T::from(1.0).unwrap()];

    if let Some(pre) = prenodes {
        s = [pre].concat()
    }

    let mut fs: Vec<T> = s.iter().map(|&si| f(si)).collect();
    let mut fx: Vec<T>;

    let mut bests: Vec<T> = Vec::with_capacity(degree as usize);
    let mut bestfs: Vec<T> = Vec::with_capacity(degree as usize);
    let mut bestw: Col<T> = Col::zeros(zero());

    let mut bestm: i32 = zero();
    let mut besterr: f64 = f64::INFINITY;

    let mut x: Vec<T>;
    let mut r: Vec<T>;
    let mut c: Mat<T>;
    let mut a: Mat<T>;

    loop {
        let m: i32 = s.len() as i32;

        x = refine(&s, std::cmp::max(refinement, 16 - m) as usize);

        fx = x.par_iter().map(|&xi| f(xi)).collect();

        c = Mat::from_fn(x.len(), s.len(), |i, j| (*x.get(i).unwrap() - *s.get(j).unwrap()).recip());
        a = Mat::from_fn(fx.len(), fs.len(), |i, j| (*fx.get(i).unwrap() - *fs.get(j).unwrap()) * *c.get(i, j));

        let w = calculate_weights(&a);

        r = Vec::with_capacity(c.nrows());
        for i in 0..c.nrows() {
            let mut weighted_sum: T = zero();
            let mut unweighted_sum: T = zero();
            for j in 0..c.ncols() {
                weighted_sum += *c.get(i, j) * *w.get(j) * *fs.get(j).unwrap();
                unweighted_sum += *c.get(i, j) * *w.get(j);
            }
            r.push(weighted_sum / unweighted_sum);
        }

        let (pos, err) = fx
            .iter()
            .zip(r.iter())
            .map(|(&fx, &r)| (fx - r).abs())
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(&b).unwrap()).unwrap();

        let err = err.to_f64().unwrap().abs();

        if (err.abs() < besterr) && err.is_finite() && err.abs() != 0.0 && m > 2 {
            bestm = m;
            besterr = err;
            bests = s.clone();
            bestfs = fs.clone();
            bestw = w.clone();
        }

        if besterr <= tol || m == degree || m - bestm == lookahead {
            break;
        }

        s.push(*x.get(pos).unwrap());
        fs.push(*fx.get(pos).unwrap());

    }

    Barycentric {
        nodes: bests.clone(),
        values: bestfs.clone(),
        weights: bestw.iter().map(|f| *f).collect::<Vec<T>>(),
        weightsxvals: bestw
            .iter()
            .zip(bestfs.iter())
            .map(|(&w, &f)| w * f )
            .collect(),
        error: besterr,
    }
}

fn refine<T>(s: &[T], p: usize) -> Vec<T>
where T: ComplexFloat + ComplexField
{
    let mut s = s.to_vec();
    s.sort_by(|a, b| a.re().partial_cmp(&b.re()).unwrap());
    let mut result = Vec::new();

    for i in 0..s.len() - 1 {
        let s_val = s[i];
        let diff_val = s[i + 1] - s_val;
        let d: Vec<T> = (1..=p).map(|j| T::from(j).unwrap() / (T::from(p + 1).unwrap())).collect();

        for d_val in d {
            result.push(s_val + d_val * diff_val);
        }
    }
    result
}