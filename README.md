<p align="center">
<img src="img/aaalligator3.jpg" />

An attempt at implementing a rust variant of the methods described in "AAA rational approximation on a continuum" 2023 (Toby Driscoll, Yuji Nakatsukasa, Lloyd N. Trefethen) 10.48550/arXiv.2305.03677 https://arxiv.org/abs/2305.03677
 
Forgive me am i fairly new to programming so i'm pretty sure this code is full of bugs and things i can't work out how to fix yet. It's also missing a fair few features mentioned in the paper like lawson iterations and domains other than [-1, 1] and there are a few logic issues in the code but it otherwise appears to mostly work! :)

I ended up having to make two variants of the code aaaxf (floating point functions) and aaaxc (complex functions) due to not knowing if there was a way to force rust to take either as an input like matlab and julia both often seem to have automatic type conversion.
</p>

## usage example
```
use aaalligator::aaaxf;

fn f(x: &f64) -> f64 {
    f64::exp(f64::cos(4.0 * x) - f64::sin(3.0 * x))
}

let a = aaaxf(&f, &150, &0, &(1000.0*f64::EPSILON)).unwrap();
let r = aaaxf(&f, &150, &0, &(1000.0*f64::EPSILON)).unwrap().r;


    println!("f(x) - r(x) {:?}", f(&0.5) - r(&0.5));

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
```
```
f(0.5) - r(0.5) 5.570544026056723e-14

    Poles:
-1.4251836349051166-0.55103515549739i
-1.4251836349051166+0.55103515549739i
-1.2793065037823754+-0i
-0.3443656711931881-0.6997077840817748i
-0.3443656711931881+0.6997077840817748i
....

Zeros:
-1.279249421913843+-0i
-0.9766906199466148-0.6647089658203892i
-0.9766906199466148+0.6647089658203892i
-0.8086181116012913-0.6725072603508072i
-0.8086181116012913+0.6725072603508072i
-0.654947258895926-0.7054914738178846i
-0.654947258895926+0.7054914738178846i
...

final error 3.643974011424689e-12
```

## bugs / todo
- lawson iterations not added at all
- "imaginary axis or right half-plane" (aaai) and "unit circle or disk" (aaaz) not implemented
- mystery very larges poles/zeros after solving eigenvalues that don't appear in julia/matlab
- sometimes results are a bit off from other versions? something to do with the loop exit conditions maybe
- after adding callable r the parallel examples broke and i don't know how to fix :\
- the last iteration inside the loop doesn't seem to get saved(?) and the final poles comes out one short? urgh
- debug builds compiled without opt-level = 1 or higher crash with stack overflow on windows (at least for me)
- ideally this package would be made closer in implementation to the cleaner/better version of the AAA at https://github.com/complexvariables/RationalFunctionApproximation.jl
- would be nice to implement a domain colouring method for visualisation of the functions and graphing functionality too
