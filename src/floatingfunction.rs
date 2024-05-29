use faer::{complex_native::c64, solvers::SpSolver, ComplexField, Mat};
use crate::AAAxfResult;

pub fn aaaxf<'a, F>(f: &'a F, degree: &'a usize, _lawson: &'a usize, tol: &'a f64) -> Option<AAAxfResult<'a>>
where
    F: Fn(&f64) -> f64,
{
    let mut s: Vec<f64> = vec![-1.0, 1.0]; // vector of support points
    let f0: Vec<f64> = s.iter().map(|&x| f(&x)).chain(xs(&s, &10).iter().map(|&x| f(&x))).collect();
    let err: f64 = f0.iter().map(|x| x.powi(2)).sum::<f64>().sqrt(); // ...or degree==0

    if err == 0.0 || *degree == 0 {
        return None;
    }

    let mut best: Option<(usize, Vec<f64>, Vec<f64>)> = None; // track the best so far
    let mut err: Vec<f64> = vec![]; // convergence statistics
    let mut nbad: Vec<usize> = vec![];

    loop { //MAIN AAA LOOP

    let m: i32 = s.len() as i32;
    let x: Vec<f64> = xs(&s, &(std::cmp::max(3, 16 - m) as usize)); // vector of test points

        let fx: Vec<f64> = x.iter().map(|&xi| f(&xi)).collect(); // evaluate f
        let fs: Vec<f64> = s.iter().map(|&si| f(&si)).collect();

        let c: Mat<f64> = Mat::from_fn(x.len(), s.len(), |i, j| 1.0 / (x[i] - s[j])); // Cauchy matrix
        let mut a: Mat<f64> = Mat::from_fn(fx.len(), fs.len(), |i, j| fx[i] - fs[j]); // Loewner matrix

        for i in 0..c.nrows() {
            for j in 0..c.ncols() {
                a.write(i, j, a.get(i, j) * c.get(i, j)) 
            }
        }

        let v: Mat<f64> = a.svd().v().to_owned(); // SVD
        let w: Vec<f64> = v.col_as_slice(v.ncols() - 1).to_vec(); // end col
        let wxfs: Vec<f64> = w.iter().zip(fs.iter()).map(|(&w, &fs)|w * fs).collect::<Vec<f64>>();

        let mut r: Vec<f64> = Vec::with_capacity(c.nrows());
        for i in 0..c.nrows() {
            let mut weighted_sum = 0.0;
            let mut unweighted_sum = 0.0;

            for j in 0..c.ncols() {
                weighted_sum += c.get(i, j) * wxfs[j];
                unweighted_sum += c.get(i, j) * w[j];
            }

            if unweighted_sum != 0.0 {
                r.push(weighted_sum / unweighted_sum);
            } else {
                r.push(0.0);
            }
        };

        let err_val: f64 = fx.iter().zip(r.iter()).map(|(fx, r)| (fx - r).abs()).max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap(); // track max error
        err.push(err_val);

        let pol: Vec<c64> = poles(&s, &w); //poles of this approximant

        let bad = pol.iter().filter(|&p| p.im == 0.0 && p.abs() <= 1.0).count(); // flag bad poles
        nbad.push(bad); // track number of bad poles

        let fmax = fs.iter().chain(fx.iter()).map(|x| x.abs()).fold(f64::NEG_INFINITY, f64::max); // set scale of f

        if best.is_none() || (bad == 0 && err.last().unwrap() < &err[best.as_ref().unwrap().0 - 1]) { 
                best = Some((m as usize, s.clone(), w));
            }


        if let Some((best_m, _, _)) = best {
            if err.len() >= best_m { // Ensure err has enough elements 
                let is_low = err[best_m - 1] / fmax < 1e-2;

                if (bad > 0 && err.last().unwrap() / fmax <= *tol)
                    || m == (degree + 1).try_into().unwrap()
                    || (m - best_m as i32 >= 10 && is_low)
                {
                    break;
                }
            }
        }

        let j = (0..fx.len()).max_by(|&i, &j| (fx[i] - r[i]).abs().partial_cmp(&(fx[j] - r[j]).abs()).unwrap()).unwrap();
        s.push(x[j]); //...and include it
    }


    let (_, s, w) = best.unwrap();

    let fs: Vec<f64> = s.iter().map(|&si| f(&si)).collect();

    let s_clone = s.clone(); // Clone s before moving it into the closure

    let (pol, res, zer) = prz(&s_clone, &fs, &w); // poles, residues, zeros

    let xx: Vec<f64> = xs(&s, &30);

    let r = move |x: &f64| {
        evaluator(&s, &fs, &w)(*x) // Call the evaluator closure
    };

    let ee: Vec<f64> = xx.iter().map(|&x| f(&x) - r(&x)).collect(); // COMPUTE ERROR AND PLOT
    let final_err = ee.iter().map(|e| e.abs()).max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap(); // Find maximum absolute error

    Some(AAAxfResult {
        r: Box::new(r),
        poles: pol,
        residues: res,
        zeros: zer,
        final_error: final_err,
        support_points: s_clone,
        vals: xx,
    })
}

fn evaluator<'a>(zj: &'a [f64], fj: &'a [f64], wj: &'a [f64]) -> Box<dyn Fn(f64) -> f64 + 'a> {
    Box::new(move |z: f64| { 
        if z.is_infinite() {
            return wj.iter().zip(fj.iter()).map(|(&wj_k, &fj_k)| wj_k * fj_k).sum::<f64>() / wj.iter().sum::<f64>();
        }

        if let Some(k) = zj.iter().position(|&zj_k| z == zj_k) {
            fj[k]
        } else {
            let c: Vec<f64> = zj.iter().map(|&zj_k| 1.0 / (z - zj_k)).collect();
            let numerator: f64 = c
                .iter()
                .zip(wj.iter())
                .zip(fj.iter())
                .map(|((&c_k, &wj_k), &fj_k)| c_k * wj_k * fj_k)
                .sum();
            let denominator: f64 = c.iter().zip(wj.iter()).map(|(&c_k, &wj_k)| c_k * wj_k).sum();
            numerator / denominator
        }
    })
}

fn prz(zj: &[f64], fj: &[f64], wj: &[f64]) -> (Vec<c64>, Vec<c64>, Vec<c64>) {
    let pol = poles(zj, wj);
    let res = residues(&pol, zj, fj, wj);
    let zer = zeros(zj, fj, wj);
    (pol, res, zer) 
}

fn poles(zj: &[f64], wj: &[f64]) -> Vec<c64>
{
    let wj = wj.iter().filter(|&f| *f != 0.0).cloned().collect::<Vec<f64>>();
    let m = wj.len();
    let mut etemp: Mat<f64> = Mat::zeros(m + 1, m + 1);

    for i in 0..m {
    etemp.write( i + 1, 0, 1.0);
    etemp.write( 0, i + 1, wj[i]);
    etemp.write(i + 1, i + 1, zj[i]);
    }

    let mut b: Mat<f64> = Mat::identity(m + 1, m + 1);
    b.write(0, 0, 0.0);

    let e: Mat<f64> = etemp.qr().solve(&b);

    let mut pol: Vec<c64> = e.eigenvalues()
    .iter()
    .map(|&f: &c64| c64::faer_one() / f)
    .filter(|&f: &c64| f.abs() < 1e4)
    .collect();

    pol.sort_by(|a, b|a.re().partial_cmp(&b.re()).unwrap());

    pol
}

fn residues(pol: &[c64], zj: &[f64], fj: &[f64], wj: &[f64]) -> Vec<c64> {

    let n = |t: c64| zj.iter().zip(fj.iter()).zip(wj.iter())
        .map(|((x, y), w)| *y * *w / (t - c64::new(*x, 0.0)))
        .sum::<c64>();

    let ddiff = |t: c64| zj.iter().zip(wj.iter())
        .map(|(&x, &w)| -w / (t - c64::new(x, 0.0)).powi(2))
        .sum::<c64>();

    let res: Vec<c64> = pol.iter()
        .map(|&p| n(p) / ddiff(p))
        //.filter(|&f: &c64| f != c64::ZERO && f.abs() < 1e4)
        .collect();

    res
}

fn zeros(zj: &[f64], fj: &[f64], wj: &[f64]) -> Vec<c64> {

    let m = wj.len();
    let mut e: Mat<f64> = Mat::zeros(m + 1, m + 1);

    for i in 0..m {
        e.write(i + 1, 0, 1.0);
        e.write(0, i + 1, wj[i] * fj[i]); 
        e.write(i + 1, i + 1, zj[i]);
    }

    let mut b: Mat<f64> = Mat::identity(m + 1, m + 1);
    b.write(0, 0, 0.0);

    let zer = e.qr().solve(&b);

    let mut zerfiltered: Vec<c64> = zer.eigenvalues()
        .iter()
        .map(|&f: &c64| c64::faer_one() / f)
        .filter(|&f: &c64| f != c64::faer_one() && f.abs() < 1e4)
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