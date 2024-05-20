pub mod floatingfunction;
pub mod complexfunction;
pub use crate::floatingfunction::aaaxf;
pub use crate::complexfunction::aaaxc;
use num_complex::Complex64;

pub struct AaaxResult {
    pub poles: Vec<Complex64>,
    pub residues: Vec<Complex64>,
    pub zeros: Vec<Complex64>,
    pub final_error: f64,
    pub support_points: Vec<f64>,
    pub vals: Vec<f64>
}