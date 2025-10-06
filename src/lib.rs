pub mod barycentric;

pub mod types;

pub mod aaa;
pub use aaa::*;

pub mod util;
pub use util::*;

mod continuum_tests;

/// struct containing the rational approximation
#[derive(Debug, Clone, Default)]
pub struct Barycentric<T> {
    pub nodes: Vec<T>,
    pub values: Vec<T>,
    pub weights: Vec<T>,
    pub weightsxvals: Vec<T>,
    pub error: f64,
}