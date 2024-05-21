use num_complex::{Complex64, ComplexFloat};
use faer::{solvers::SpSolver, Mat};
use crate::AAAxcResult;

pub fn aaaxc<'a, F>(f: &'a F, degree: &'a usize, _lawson: &'a usize, tol: &'a f64) -> Option<AAAxcResult<'a>>
where
    F: Fn(&Complex64) -> Complex64,
{
    let mut s: Vec<f64> = vec![-1.0, 1.0]; // vector of support points
    let f0: Vec<Complex64> = s.iter().map(|&x| f(&x.into())).chain(xsc(&s, &10).iter().map(|&x| f(&x))).collect(); //# check for constant function...
    let err: f64 = f0.iter().map(|x| x.abs().powi(2)).sum::<f64>().sqrt(); // ...or degree==0

    if err == 0.0 || *degree == 0 {
        return None;
    }

    let mut best: Option<(usize, Vec<f64>, Vec<Complex64>)> = None; // track the best so far
    let mut err: Vec<f64> = vec![]; // convergence statistics
    let mut nbad: Vec<usize> = vec![];

    loop { //MAIN AAA LOOP
        let m: i32 = s.len() as i32;

        let x: Vec<f64> = xs(&s, &(std::cmp::max(3, 16 - m) as usize)); // vector of test points

        let fx: Vec<Complex64> = x.iter().map(|&xi| f(&xi.into())).collect(); // evaluate f
        let fs: Vec<Complex64> = s.iter().map(|&si| f(&si.into())).collect();

        let c: Mat<f64> = Mat::from_fn(x.len(), s.len(), |i, j| 1.0 / (x[i] - s[j])); // Cauchy matrix
        let mut a: Mat<Complex64> = Mat::from_fn(fx.len(), fs.len(), |i, j| fx[i] - fs[j]); // Loewner matrix

        for i in 0..a.nrows() {
            for j in 0..a.ncols() {
                a.write(i, j, a.read(i, j) * c.read(i, j)) 
            }
        }

        let v: Mat<Complex64> = a.svd().v().to_owned(); // SVD
        let w: Vec<Complex64> = v.col_as_slice(v.ncols() - 1).re.iter().zip(v.col_as_slice(v.ncols() - 1).im).map(|(re, im)|Complex64::new(*re, *im)).collect::<Vec<Complex64>>().to_vec(); // right col
        let wxfs: Vec<Complex64> = w.iter().zip(fs.iter()).map(|(&w, &fs)|w * fs).collect::<Vec<Complex64>>();

        let mut r: Vec<Complex64> = Vec::with_capacity(c.nrows());
        for i in 0..c.nrows() {
            let mut weighted_sum: Complex64 = Complex64::ZERO;
            let mut unweighted_sum: Complex64 = Complex64::ZERO;

            for j in 0..c.ncols() {
                weighted_sum += c.read(i, j) * wxfs[j];
                unweighted_sum += c.read(i, j) * w[j];
            }

            if unweighted_sum != Complex64::ZERO {
                r.push(weighted_sum / unweighted_sum);
            } else {
                r.push(Complex64::ZERO);
            }
        };

        let err_val: f64 = fx.iter().zip(r.iter()).map(|(fx, r)| (fx - r).abs()).max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap(); // track max error
        err.push(err_val);

        let pol: Vec<Complex64> = poles(&s, &w); //poles of this approximant

        let bad = pol.iter().filter(|&p| p.im == 0.0 && p.abs() <= 1.0).count(); // flag bad poles
        nbad.push(bad); // track number of bad poles

        let fmax = fs.iter().chain(fx.iter()).map(|x| x.abs()).fold(f64::NEG_INFINITY, f64::max); // set scale of f

        if best.is_none() || (bad == 0 && err.last().unwrap() < &err[best.as_ref().unwrap().0 - 1]) { 
                best = Some((m as usize, s.clone(), w));
            }

        if let Some((best_m, _, _)) = best {
            if err.len() >= best_m { // Ensure err has enough elements 
                let is_low = err[best_m - 1] / fmax < 1e-2;

                if (bad == 0 && err.last().unwrap() / fmax <= *tol)
                    || m == (degree + 1).try_into().unwrap()
                    || (m - best_m as i32 >= 10 && is_low)
                {
                    break;
                }
            }
        }

        // find next support point...
        let j = (0..x.len()).max_by(|&i, &j| (fx[i] - r[i]).abs().partial_cmp(&(fx[j] - r[j]).abs()).unwrap()).unwrap();
        s.push(x[j]); //...and include it
    }

    let (_, s, w) = best.unwrap();

    let fs: Vec<Complex64> = s.iter().map(|&si| f(&si.into())).collect();

    let s_clone = s.clone(); // Clone s before moving it into the closure

    let (pol, res, zer) = prz(&s, &fs, &w); // poles, residues, zeros

    let xx: Vec<f64> = xs(&s, &30);

    let r = move |x: &Complex64| {
        evaluator(&s, &fs, &w)(*x) // Call the evaluator closure
    };

    let ee: Vec<Complex64> = xx.iter().map(|&x| f(&x.into()) - r(&x.into())).collect(); // COMPUTE ERROR AND PLOT
    let final_err = ee.iter().map(|e| e.abs()).max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap(); // Find maximum absolute error

    Some(AAAxcResult {
        r: Box::new(r),
        poles: pol,
        residues: res,
        zeros: zer,
        final_error: final_err,
        support_points: s_clone,
        vals: xx,
    })
}

fn evaluator<'a>(zj: &'a [f64], fj: &'a [Complex64], wj: &'a [Complex64]) -> Box<dyn Fn(Complex64) -> Complex64 + 'a> {
    Box::new(move |z: Complex64| { 
        if z.is_infinite() {
            return wj.iter().zip(fj.iter()).map(|(&wj_k, &fj_k)| wj_k * fj_k).sum::<Complex64>() / wj.iter().sum::<Complex64>();
        }

        if let Some(k) = zj.iter().position(|&zj_k| z == zj_k.into()) {
            fj[k]
        } else {
            let c: Vec<Complex64> = zj.iter().map(|&zj_k| 1.0 / (z - zj_k)).collect();
            let numerator: Complex64 = c
                .iter()
                .zip(wj.iter())
                .zip(fj.iter())
                .map(|((&c_k, &wj_k), &fj_k)| c_k * wj_k * fj_k)
                .sum();
            let denominator: Complex64 = c.iter().zip(wj.iter()).map(|(&c_k, &wj_k)| c_k * wj_k).sum();
            numerator / denominator
        }
    })
}

fn prz(zj: &[f64], fj: &[Complex64], wj: &[Complex64]) -> (Vec<Complex64>, Vec<Complex64>, Vec<Complex64>) {
    let pol = poles(zj, wj);
    let res = residues(&pol, zj, fj, wj);
    let zer = zeros(zj, fj, wj);
    (pol, res, zer) 
}

fn poles(zj: &[f64], wj: &[Complex64]) -> Vec<Complex64> {
    let wj = wj.iter().filter(|&f| *f != Complex64::ZERO).cloned().collect::<Vec<Complex64>>();
    let m = wj.len();
    let mut etemp: Mat<Complex64> = Mat::zeros(m + 1, m + 1);

    for i in 0..m {
        etemp.write(i + 1, 0, Complex64::ONE);
        etemp.write(0, i + 1, wj[i]);
        etemp.write(i + 1, i + 1, zj[i].into());
    }

    let mut b: Mat<Complex64> = Mat::identity(m + 1, m + 1);
    b.write(0, 0, Complex64::ZERO);

    let e: Mat<Complex64> = etemp.qr().solve(&b);

    let mut pol: Vec<Complex64> = e.eigenvalues()
        .iter()
        .map(|&f: &Complex64| Complex64::ONE / f)
        .filter(|&f: &Complex64| f.abs() < 1e4)
        .collect();

    pol.sort_by(|a, b| a.re().partial_cmp(&b.re()).unwrap());

    pol
}

fn residues(pol: &[Complex64], zj: &[f64], fj: &[Complex64], wj: &[Complex64]) -> Vec<Complex64> {

    let n = |t: Complex64| zj.iter().zip(fj.iter()).zip(wj.iter())
        .map(|((x, y), w)| *y * *w / (t - Complex64::new(*x, 0.0)))
        .sum::<Complex64>();

    let ddiff = |t: Complex64| zj.iter().zip(wj.iter())
        .map(|(&x, &w)| -w / (t - Complex64::new(x, 0.0)).powi(2))
        .sum::<Complex64>();

    let res: Vec<Complex64> = pol.iter()
        .map(|&p| n(p) / ddiff(p))
        //.filter(|&f: &Complex64| f != Complex64::ZERO && f.abs() < 1e14)
        .collect();

    res
}

fn zeros(zj: &[f64], fj: &[Complex64], wj: &[Complex64]) -> Vec<Complex64> {

    let m = wj.len();
    let mut e: Mat<Complex64> = Mat::zeros(m + 1, m + 1);

    for i in 0..m {
        e.write(i + 1, 0, Complex64::ONE);
        e.write(0, i + 1, wj[i] * fj[i]); 
        e.write(i + 1, i + 1, zj[i].into());
    }

    let mut b: Mat<Complex64> = Mat::identity(m + 1, m + 1);
    b.write(0, 0, Complex64::ZERO);

    let zer = e.qr().solve(&b);

    let mut zerfiltered: Vec<Complex64> = zer.eigenvalues()
        .iter()
        .map(|&f: &Complex64| Complex64::ONE / f)
        .filter(|&f: &Complex64| f != Complex64::ZERO && f.abs() < 1e4)
        .collect();

    zerfiltered.sort_by(|a, b| a.re().partial_cmp(&b.re()).unwrap());

    zerfiltered
}

fn xs(s: &[f64], p: &usize) -> Vec<f64> {
    let mut s = s.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut result = Vec::new(); 
    for i in 0..s.len() - 1 {
        let s_val = s[i];
        let diff_val = s[i + 1] - s_val;

        let d: Vec<f64> = (1..=*p).map(|j| j as f64 / ((p + 1) as f64)).collect();

        for d_val in d {
            result.push(s_val + d_val * diff_val);
        }
    }

    result
}

fn xsc(s: &[f64], p: &usize) -> Vec<Complex64> {
    let s: &mut [f64] = &mut s.to_owned();
    s.sort_by(|a, b| a.abs().clone().partial_cmp(&b.abs().clone()).unwrap());

    let d: Vec<f64> = (1..=*p).map(|i| i as f64 / (p + 1) as f64).collect();
    let mut result: Vec<Complex64> = Vec::new();
    for i in 0..(s.len() - 1) {
        for &j in d.iter() {
            result.push((s[i] + j * (s[i + 1] - s[i])).into());
        }
    }
    result
}