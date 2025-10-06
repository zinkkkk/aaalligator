use std::f64::consts::TAU;

use aaalligator::*;
use num_complex::{Complex64, ComplexFloat};

fn main() {

    let r = f_f64.into_bary();

    r.print_info();
    r.przr().print_przr_info();
    //r.draw();

    // let r = kite.into_bary();
    // r.print_info();
    // r.przr().print_przr_info();
    // r.draw_with_pz();

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