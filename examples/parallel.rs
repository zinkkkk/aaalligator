use aaalligator::{aaaxc, aaaxf, AaaxResult};
use num_complex::Complex64;
use rayon::prelude::*;
use std::sync::mpsc::channel;

fn main() {

    // run 100x

    println!("f = x -> exp(cos(4.0 + offset * x) - sin(3.0 * x))");
    let (sender, receiver) = channel();
    for i in 0..100 {
        let sender = sender.clone();
        std::thread::spawn(move || {
            let offset = i as f64 * 0.01;
            let f_with_offset = |x: &f64| f(x, offset); 
            let result = aaaxf(&f_with_offset, &150, &0, &(1000.0 * f64::EPSILON)); 
            sender.send(result).unwrap(); 
        });
    }

    let results: Vec<AaaxResult> = receiver.iter().take(100).collect();

    for result in results {
        println!("{:?}", result.poles);
        println!("{:?}", result.zeros);
        println!("{:?}", result.final_error);
    }
    

    // with rayon

    println!("fc = (z, offset) -> (z^3 - 1) / sin(z - Complex(0.9 + offset, 1.0))");
    let rayon: Vec<AaaxResult> = (0..100)
        .into_par_iter() 
        .map(|i| {
            let offset = i as f64 * 0.01;
            aaaxc(&|z| fc(z, &offset), &150, &0, &(1000.0 * f64::EPSILON))
        })
        .collect();

    for result in rayon {
        println!("vals len: {:?}", result.poles);
        println!("vals: {:?}", result.zeros);
        println!("vals: {:?}", result.final_error);
    }

}

fn f(x: &f64, offset: f64) -> f64 {
    f64::exp(f64::cos(4.0 + offset * x) - f64::sin(3.0 * x))
}

fn fc(z: &Complex64, offset: &f64) -> Complex64 {
    (z.powi(3) - Complex64::ONE) / Complex64::sin(z - (Complex64::new(0.9 + offset, 1.0)))
}

