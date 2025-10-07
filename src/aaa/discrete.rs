use crate::*;

use faer::Col;
use faer::prelude::*;
use faer::traits::ComplexField;
use num_complex::ComplexFloat;
use num_traits::NumAssignOps;
use num_traits::NumOps;
use num_traits::ToPrimitive;

/// AAA approximation from sets of data contrained in Vec type
pub fn aaa_discreet<T>(
    f: Vec<T>,
    z: Vec<T>,
    degree: usize,
    tol: f64,
    lookahead: usize
) -> Barycentric<T>
where T: ComplexFloat + ComplexField + NumOps<T> + NumAssignOps<T>, f64: From<<T as ComplexField>::Real>
{
    let mut m: usize;
    let fzlen: usize = z.len();

    let mut zj = Vec::with_capacity(degree);
    let mut fj = Vec::with_capacity(degree);
    let mut wj;

    let mut bestzj = Vec::with_capacity(degree);
    let mut bestfj = Vec::with_capacity(degree);
    let mut bestwj = Col::zeros(degree);

    let mut bestm = 0;
    let mut besterr = f64::INFINITY;
    let mut errs: Vec<f64>;
    let mut err: f64;

    let mut c: Mat<T>;
    let mut a: Mat<T>;

    let mut rr = Vec::with_capacity(fzlen);
    let mut j: Vec<usize> = (0..fzlen).collect();
    let mut jj: usize;
    let mut r = vec![ColRef::from_slice(f.as_slice()).sum() / T::from(fzlen).unwrap(); fzlen];

    loop {

        rr.clear();
        for i in &j {
            rr.push((*f.get(*i).unwrap() - *r.get(*i).unwrap()).abs());
        }

        jj = rr
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).unwrap())
            .unwrap().0;

        zj.push(*z.get(*j.get(jj).unwrap()).unwrap());
        fj.push(*f.get(*j.get(jj).unwrap()).unwrap());
        j.remove(jj);
        m = zj.len();

        c = Mat::from_fn(z.len(), zj.len(), |i, j| (*z.get(i).unwrap() - *zj.get(j).unwrap()).recip());
        a = Mat::from_fn(f.len(), fj.len(), |i, j| (*f.get(i).unwrap() - *fj.get(j).unwrap()) * *c.get(i, j));

        let mut a2: Mat<T> = Mat::zeros(j.len(), zj.len());
        for i in 0..m {
            for (k, j) in j.iter().enumerate() {
                *a2.get_mut(k, i) = *a.get(*j, i);
            }
        }

        wj = calculate_weights(&a2);

        let mut wjfj: Col<T> = Col::zeros(wj.nrows());
        wjfj.iter_mut().zip(wj.iter().zip(fj.iter())).for_each(|(a, (i, j))| *a = *i * *j);

        let n: Vec<T> = (&c * wjfj).to_owned().try_as_col_major().unwrap().as_slice().to_vec();
        let d: Vec<T> = (&c * &wj).to_owned().try_as_col_major().unwrap().as_slice().to_vec();

        r = n
            .as_slice()
            .iter()
            .zip(d.as_slice().iter())
            .map(|(n_val, d_val)| *n_val / *d_val )
            .collect();

        errs = f
            .as_slice()
            .iter()
            .zip(r.as_slice().iter())
            .map(|(ff, rr)| ((ff.abs() - rr.abs()).to_f64().unwrap()).abs())
            .collect();

        errs = errs.iter().cloned().filter(|&a| a.is_normal()).collect();
        err = RowRef::from_slice(errs.iter().as_slice()).norm_max().into();

        if err < besterr {
            bestm = m as i32;
            besterr = err;
            bestzj = zj.to_owned();
            bestfj = fj.to_owned();
            bestwj = wj.clone();
        }

        if besterr <= tol || m == degree || (bestm as usize) + (lookahead as usize) == m || m + 1 == fzlen
        { break; }
    }

    Barycentric {
        nodes: bestzj.clone(),
        values: bestfj.clone(),
        weights: bestwj.clone().try_as_col_major().unwrap().as_slice().to_vec(),
        weightsxvals: bestwj
            .iter()
            .zip(bestfj.iter())
            .map(|(w, f)| *w * *f)
            .collect(),
        error: besterr,
    }
}
