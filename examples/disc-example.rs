use std::f64::consts::TAU;

use aaalligator::*;
use iter_num_tools::lin_space;
use num_complex::{Complex64, ComplexFloat};

fn main() {

    let z = lin_space(-1.0..=1.0, 500).collect::<Vec<f64>>();
    let f = z.iter().map(|&i| f_f64(i)).collect::<Vec<f64>>();

    f.clone().interpolate_n(1000).draw();

    let r = (f, z).into_bary();
    r.draw();

    let z = lin_space(Complex64::from(-1.0)..=Complex64::from(1.0), 500).collect::<Vec<Complex64>>();
    let f = z.iter().map(|&i| kite(i)).collect::<Vec<Complex64>>();

    let r = f.into_bary();
    r.draw_with_pz();

}

fn f_f64(x: f64) -> f64 {
    ((4.0 * x).cos() - (3.0 * x).sin()).exp()
}

fn kite(t: Complex64) -> Complex64 {
    let input_param = t.re();

    let t = (input_param + 1.0) * 0.5 * TAU;

    let x = 2.2 * t.cos() + 1.25 * (2.0 * t).cos() - 1.25;
    let y = 3.0 * t.sin();

    Complex64::new(x, y)
}