#[cfg(test)]
mod tests {
    use crate::*;
    use crate::util::consts_internal::*;

    #[test]
    fn basic_functions() {

        let f = |x: f64| (1.05 - x).sin();
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 2e-13);

        let f = |x: f64| (-1.0 / (x * x)).exp();
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        println!("dif {:?}", dif / (NEGONE_TO_ONE.len() as f64));
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 4e-13);

        let f = |x: f64| (-100.0 * x * x).exp();
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 2e-13);

        let f = |x: f64| (-10.0 / (1.2 - x)).exp();
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 1e-12);

        let f = |x: f64| 1.0 / (1.0 + (100.0 * (x + 0.5)).exp());
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 2e-13);

        // i dunno whats the go with this one lol;
        // let f = |x: f64| (100.0 * x).sin() * (-10.0 * x * x).exp();
        // let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        // let r = approx.evaluate_handle();
        // let mut dif = 0.0;
        // for i in NEGONE_TO_ONE {
        //     dif += (r(i) - f(i)).abs();
        // }
        // println!("diff {}", dif);
        // println!("deg {}", approx.nodes.len());
        // assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 1e-11);

        let f = |x: f64| x.abs();
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 1e-8);

        let f = |x: f64| (x - 0.95).abs();
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 1e-6);
    }

    #[test]
    fn low_accuracy() {
        let f = |x: f64| (3.0 * x).exp();
        let approx = aaa_continuum(f, None, 150, 0, 1e-4, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) > 1e-8);
    }

    #[test]
    fn vertical_scaling() {
        let f = |x: f64| 1e100 * x.sin();
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += r(i) - f(i);
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 2e-13);

        println!("start ver");
        let f = |x: f64| 1e-100 * x.cos();
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        println!("diff {}", dif);
        println!("deg {}", approx.nodes.len());
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 2e-13);
    }

    #[test]
    fn polynomials_and_reciprocals() {
        let f = |x: f64| x;
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 2e-13);

        let f = |x: f64| x + x * x;
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 2e-13);

        let f = |x: f64| x + x * x * x;
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 2e-13);

        let f = |x: f64| 1.0 / (1.1 + x);
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 2e-13);

        let f = |x: f64| 1.0 / (3.0 + x + x * x);
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 2e-13);

        let f = |x: f64| 1.0 / (1.01 + x * x * x);
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 2e-13);
    }

    #[test]
    fn specified() {
        let f = |x: f64| x;
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 2e-13);

        let f = |x: f64| x + x * x;
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 2e-13);

        let f = |x: f64| x + x * x * x;
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 2e-13);

        let f = |x: f64| 1.0 / (1.1 + x);
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 2e-13);

        let f = |x: f64| 1.0 / (3.0 + x + x * x);
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 2e-13);

        let f = |x: f64| 1.0 / (1.01 + x * x * x);
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 2e-13);

        let f = |x: f64| (100.0 * x).tanh();
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += (r(i) - f(i)).abs();
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 2e-13);

        let f = |x: f64| (100.0 * (x - 0.2)).tanh();
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += r(i) - f(i);
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 2e-13);

        let f = |x: f64| x.exp();
        let approx = aaa_continuum(f, None, 150, 0, 1000.0 * f64::EPSILON, 10, 3);
        let r = approx.evaluate_handle();
        let mut dif = 0.0;
        for i in NEGONE_TO_ONE {
            dif += r(i) - f(i);
        }
        assert!(dif / (NEGONE_TO_ONE.len() as f64) <= 2e-13);
    }
}