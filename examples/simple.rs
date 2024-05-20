use aaalligator::{aaaxc, aaaxf};
use num_complex::Complex64;

fn main() {

    // floating point function approximations

    println!("f = x -> exp(cos(4x) - sin(3x))");
    let a = aaaxf(&f, &150, &0, &(1000.0*f64::EPSILON));
    println!("Poles:");
    for i in a.poles {
    println!("{}", i)    
    }
    println!();
    println!("Zeros:");
    for i in a.zeros {
    println!("{}", i)    
    }
    println!("\nfinal error {:?}\n", a.final_error);

    println!("f = x -> tanh(10*(x - 0.1)^2)");
    let b = aaaxf(&f2, &150, &0, &(1000.0*f64::EPSILON));
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
    let c = aaaxc(&fc, &150, &0, &(1000.0*f64::EPSILON));
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
    let d = aaaxc(&fc2, &150, &0, &(1000.0*f64::EPSILON));
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

fn fc(z: &Complex64) -> Complex64 {
    (z.powi(3) - Complex64::ONE) / Complex64::sin(z - (Complex64::new(0.9, 1.0)))
}

fn fc2(z: &Complex64) -> Complex64 {
  (1.0 - (1.0 / z.powf(2.0) + z.powf(3.0))).sqrt()
}
