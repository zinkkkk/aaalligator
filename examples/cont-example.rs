use aaalligator::*;
use num_complex::{Complex64, ComplexFloat};

fn main() {

    let r = f_f64.into_bary();
    r.print_info();
    r.przr().print_przr_info();
    r.draw();

    let r = kite.into_bary();
    r.print_info();
    r.przr().print_przr_info();
    r.draw();
    r.draw_with_pz();

}

fn f_f64(x: f64) -> f64 {
    ((4.0 * x).cos() - (3.0 * x).sin()).exp()
}

fn kite(t: Complex64) -> Complex64 {
    let t = (t.re() + 1.0) * 0.5 * std::f64::consts::TAU;

    let x = 2.2 * t.cos() + 1.25 * (2.0 * t).cos() - 1.25;
    let y = 3.0 * t.sin();

    Complex64::new(x, y)
}

fn asteroid(t: Complex64) -> Complex64 {
    let t = t * std::f64::consts::PI;
    Complex64::new(t.re().cos().powi(3), t.re().sin().powi(3))
}

fn bicorn(t: Complex64) -> Complex64 {
    let t = t * std::f64::consts::PI;
    let (sin_t, cos_t) = t.re().sin_cos();

    let x = sin_t;
    
    let cos_t_sq = cos_t.powi(2);
    let sin_t_sq = sin_t.powi(2);

    let y_numerator = cos_t_sq * (2.0 + cos_t);
    let y_denominator = 3.0 + sin_t_sq;
    let y = y_numerator / y_denominator;

    Complex64::new(x, y)
}