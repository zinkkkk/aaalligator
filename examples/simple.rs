use aaalligator::{aaaxf, aaaxc};
use faer::{complex_native::c64, ComplexField};

fn main() {

    // floating point function approximations

    println!("f = x -> exp(cos(4x) - sin(3x))");
    let a = aaaxf(&f, &150, &0, &(1000.0*f64::EPSILON)).unwrap();
    let r = aaaxf(&f, &150, &0, &(1000.0*f64::EPSILON)).unwrap().r;

    println!("f: {:?}", f(&0.5));
    println!("r: {:?}", (r)(&0.5));

    println!("f(x) - r(x) {:?}", f(&0.5) - r(&0.5));

    let mut sum_difference = 0.0;
    let step: f64 = 0.01;
    let mut x = -1.0; 
    while x <= 1.0 { 
        sum_difference += f(&x) - r(&x);
        x += step;
    }
    println!("Sum of f(x) - r(x) over [-1, 1]: {:?}", sum_difference);


    println!("Poles:");
    println!("Poles len : {}", a.poles.len());
    for i in a.poles {
    println!("{}", i)    
    }
    println!();
    println!("Zeros:");
    println!("Poles len : {}", a.zeros.len());
    for i in a.zeros {
    println!("{}", i)    
    }
    println!("\nfinal error {:?}\n", a.final_error);


    println!("f = x -> tanh(10*(x - 0.1)^2)");
    let b = aaaxf(&f2, &150, &0, &(1000.0*f64::EPSILON)).unwrap();
    println!("Poles:");
    for i in b.poles {
    println!("{}", i)    
    }
    println!();
    println!("Zeros:");
    for i in b.zeros {
    println!("{}", i)    
    }
    println!("\nfinal error {:?}\n", b.final_error);

    // complex function approximations

    println!("f = z -> (z^3 - 1) / sin(z - 0.9 - 1im)");
    let c = aaaxc(&fc, &150, &0, &(1000.0*f64::EPSILON)).unwrap();
    println!("Poles:");    
    for i in c.poles {
    println!("{}", i)    
    }
    println!();
    println!("Zeros:");
    for i in c.zeros {
    println!("{}", i)    
    }
    println!("\nfinal error {:?}\n", c.final_error);

    println!("f = z -> sqrt(1.0 - (1.0 / z^2 + z^3))");
    let d = aaaxc(&fc2, &150, &0, &(1000.0*f64::EPSILON)).unwrap();
    println!("Poles:");    
    for i in d.poles {
    println!("{}", i)    
    }
    println!();
    println!("Zeros:");
    for i in d.zeros {
    println!("{}", i)    
    }
    println!("\nfinal error {:?}", d.final_error);
}

fn f(x: &f64) -> f64 {
    f64::exp(f64::cos(4.0 * x) - f64::sin(3.0 * x))
}

fn f2(x: &f64) -> f64 {
    (10.0 * (x - 0.1).powi(2)).tanh()
}

fn fc(z: &c64) -> c64 {
    (z.powi(3) - c64::faer_one()) / c64::sin(*z - (c64::new(0.9, 1.0)))
}

fn fc2(z: &c64) -> c64 {
  (1.0 - (1.0 / z.powf(2.0) + z.powf(3.0))).sqrt()
}
