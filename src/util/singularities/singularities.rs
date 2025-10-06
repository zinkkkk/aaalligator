use faer::traits::{RealField};
use num_complex::{Complex, Complex64, ComplexFloat};
use num_traits::{Float, One, Zero, float::FloatCore};
use crate::*;

/// struct containing the poles, poles residues, zeros, zeros residues of an approximation
#[derive(Debug, Clone, Default)]
pub struct Singularities<T>
{
    pub poles: Vec<Complex<T>>,
    pub residues: Vec<Complex<T>>,
    pub zeros: Vec<Complex<T>>,
}

pub trait FindSingularities<T>
{
    fn przr(&self) -> Singularities<T>;
}

impl<T> FindSingularities<T> for Barycentric<T>
where T: RealField + Float + FloatCore + ComplexFloat
{
    /// returns a struct containing poles, poles residues, zeros, zeros residues
    fn przr(&self) -> Singularities<T> {

        let poles = self.poles();
        let residues = self.residues(&poles);
        let zeros = self.zeros();

        Singularities {
            poles,
            residues,
            zeros,
        }
    }
}

impl<T: std::fmt::Debug> Singularities<T> {
    pub fn print_przr_info(&self) {
        println!("Poles");
        for i in &self.poles {
            println!("{:?}", i)
        };
        println!();

        println!("Residues");
        for i in &self.residues {
            println!("{:?}", i)
        };
        println!();

        println!("Zeros");
        for i in &self.zeros {
            println!("{:?}", i)
        };
        println!();
    }
}

impl Singularities<f64> {

    /// prints a tuple of usizes containing the number of entries in each array
    pub fn len_sizes(&self) -> (usize, usize, usize) {
        (self.poles.len(), self.residues.len(), self.zeros.len())
    }

    /// prints a string containing the numer of entries of each singularity
    pub fn len_text(&self) -> String {
        format!(
            "Poles: {}, Poles residues: {}, Zeros: {}",
            self.poles.len(),
            self.residues.len(),
            self.zeros.len(),
        )
    }

    /// sorts values by absolute. singularities by default in order matching the relitive posistion of the approximation nodes
    /// sorted by smallest to largest real part
    pub fn sort(&mut self) -> Self {

        let len: usize = self.poles.len();
        let mut indices: Vec<usize> = (0..len).collect();
        indices.sort_by(|&i, &j| self.poles[i].re().total_cmp(&self.poles[j].re()));
        let mut sorted_poles: Vec<Complex64> = Vec::with_capacity(len);
        let mut sorted_residues: Vec<Complex64> = Vec::with_capacity(len);

        for i in indices {
            sorted_poles.push(self.poles[i]);
            sorted_residues.push(self.residues[i]);
        }


        // zeros
        let len: usize = self.zeros.len();
        let mut indices: Vec<usize> = (0..len).collect();
        indices.sort_by(|&i, &j| self.zeros[i].re().total_cmp(&self.zeros[j].re()));
        let mut sorted_zeros: Vec<Complex64> = Vec::with_capacity(len);

        for i in indices {
            sorted_zeros.push(self.zeros[i]);
        }

        Self {
            poles: sorted_poles,
            residues: sorted_residues,
            zeros: sorted_zeros,
        }

    }

    /// filters out zero values and anything greater 1e-13 in size
    pub fn filter(&self) -> Self {

        let poles: Vec<Complex64> = self.poles
            .iter().copied()
            .filter(|&f: &Complex64| f != Complex64::zero() && f.abs() > 1e-13)
            .collect::<Vec<Complex64>>();

        let residues: Vec<Complex64> = self.residues
            .iter()
            .copied()
            .filter(|&f: &Complex64| f != Complex64::zero() && f.abs() > 1e-13)
            .collect::<Vec<Complex64>>();

        let zeros: Vec<Complex64> = self.zeros
            .iter().copied()
            .filter(|&f: &Complex64| f != Complex64::one() && f.abs() > 1e-13)
            .collect::<Vec<Complex64>>();

        Singularities {
            poles,
            residues,
            zeros,
        }
    }

    /// filters out zero values and anything greater than specified tolerance
    pub fn filterwithtol(&self, tol: f64) -> Self {
        let poles: Vec<Complex64> = self.poles
            .iter().copied()
            .filter(|&f: &Complex64| f != Complex64::zero() && f.abs() > tol)
            .collect::<Vec<Complex64>>();

        let residues: Vec<Complex64> = self.residues
            .iter()
            .copied()
            .filter(|&f: &Complex64| f != Complex64::zero() && f.abs() > tol)
            .collect::<Vec<Complex64>>();

        let zeros: Vec<Complex64> = self.zeros
            .iter().copied()
            .filter(|&f: &Complex64| f != Complex64::one() && f.abs() > tol)
            .collect::<Vec<Complex64>>();

        Singularities {
            poles,
            residues,
            zeros,
        }
    }
}
