pub mod floatingfunction;
pub mod complexfunction;
pub use crate::floatingfunction::aaaxf;
pub use crate::complexfunction::aaaxc;
use faer::complex_native::c64;

pub struct AAAxfResult<'a> {
    pub r: Box<dyn Fn(&f64) -> f64 + 'a>,
    pub poles: Vec<c64>,
    pub residues: Vec<c64>,
    pub zeros: Vec<c64>,
    pub final_error: f64,
    pub support_points: Vec<f64>,
    pub vals: Vec<f64>
}

pub struct AAAxcResult<'a> {
    pub r: Box<dyn Fn(&c64) -> c64 + 'a + Send>,
    pub poles: Vec<c64>,
    pub residues: Vec<c64>,
    pub zeros: Vec<c64>,
    pub final_error: f64,
    pub support_points: Vec<f64>,
    pub vals: Vec<f64>
}
