use faer::{Mat, diag::Diag, matrix_free::LinOp, unzip, zip};
use num_complex::{Complex64, ComplexFloat};
use num_traits::{One, Zero};
use crate::*;

// want to get rid of these but genertics hurt my brain legacy half finished code lol
impl Barycentric<Complex64> {

    /// helper function to sort a rational approximation by increaseing order of node value.
    /// this does not change the result when evaluated
    pub fn sort_bary(&mut self) -> Self {

        let len: usize = self.nodes.len();
        let mut indices: Vec<usize> = (0..len).collect();
        indices.sort_by(|&i, &j| self.nodes[i].re().total_cmp(&self.nodes[j].re()));

        let mut sorted_node: Vec<Complex64> = Vec::with_capacity(len);
        let mut sorted_values: Vec<Complex64> = Vec::with_capacity(len);
        let mut sorted_weights: Vec<Complex64> = Vec::with_capacity(len);
        let mut sorted_weightsxvals: Vec<Complex64> = Vec::with_capacity(len);

        for i in indices {
            sorted_node.push(self.nodes[i]);
            sorted_values.push(self.values[i]);
            sorted_weights.push(self.weights[i]);
            sorted_weightsxvals.push(self.weightsxvals[i]);
        }

        Self {
            nodes: sorted_node,
            values: sorted_values,
            weights: sorted_weights,
            weightsxvals: sorted_weightsxvals,
            error: self.error,
        }
    }

    /// returns a struct containing poles, poles residues, zeros, zeros residues
    pub fn przr(&self) -> Singularities<f64> {

        let poles: Vec<Complex64> = self.poles();
        let residues: Vec<Complex64> = self.residues(&poles);
        let zeros: Vec<Complex64> = self.zeros();

        Singularities {
            poles,
            residues,
            zeros,
        }.sort()
    }

    /// finds poles from an approximation
    pub fn poles(&self) -> Vec<Complex64> {

        let zj: &Vec<Complex64> = &self.nodes;
        let wj: &Vec<Complex64> = &self.weights;

        let wj: Vec<Complex64> = wj
            .iter()
            .filter(|&f| f.abs() != 0.0)
            .cloned()
            .collect::<Vec<Complex64>>();

        let m: usize = wj.len();

        let mut etemp: Mat<Complex64> = Mat::zeros(m + 1, m + 1);

        for i in 0..m {
            etemp[(i + 1, 0)] = Complex64::one();
            etemp[(0, i + 1)] = wj[i];
            etemp[(i + 1, i + 1)] = zj[i];
        }

        let mut b: Mat<Complex64> = Mat::identity(m + 1, m + 1);
        b[(0, 0)] = Complex64::zero();

        let gevd = etemp.generalized_eigen(b).unwrap();

        let s_a = gevd.S_a();
        let s_b = gevd.S_b();

        let mut eigs = Diag::zeros(s_a.nrows());
        zip!(&mut eigs, &s_a, &s_b).for_each(|unzip!(a, b, c)| *a = *b / *c);

        let pol: Vec<num_complex::Complex<f64>> = eigs.column_vector().try_as_col_major().unwrap().as_slice().iter().cloned().collect();

        pol
    }

    /// finds zeros from an approximation
    pub fn zeros(&self) -> Vec<Complex64> {

        let zj: &Vec<Complex64> = &self.nodes;
        let wj: &Vec<Complex64> = &self.weights;
        let wfj: &Vec<Complex64> = &self.weightsxvals;

        let m: usize = wj.len();
        let mut etemp: Mat<Complex64> = Mat::zeros(m + 1, m + 1);

        for i in 0..m {
            etemp[(i + 1, 0)] = Complex64::one();
            etemp[(0, i + 1)] = *wfj.get(i).unwrap();
            etemp[(i + 1, i + 1)] = *zj.get(i).unwrap();
        }

        let mut b: Mat<Complex64> = Mat::identity(m + 1, m + 1);
        b[(0, 0)] = Complex64::zero();

        let gevd = etemp.generalized_eigen(b).unwrap();

        let s_a = gevd.S_a();
        let s_b = gevd.S_b();

        let mut eigs = Diag::zeros(s_a.nrows());
        zip!(&mut eigs, &s_a, &s_b).for_each(|unzip!(a, b, c)| *a = *b / *c);

        let zer: Vec<num_complex::Complex<f64>> = eigs.column_vector().try_as_col_major().unwrap().as_slice().iter().cloned().collect();


        zer
    }

    /// finds residues of a singularity type
    pub fn residues(&self, polzer: &[Complex64]) -> Vec<Complex64> {

        let zj: &Vec<Complex64> = &self.nodes;
        let fj: &Vec<Complex64> = &self.values;
        let wj: &Vec<Complex64> = &self.weights;

        let n = |t: Complex64|
            zj
                .iter()
                .zip(fj.iter())
                .zip(wj.iter())
                .map(|((x, y), w)| (*y * *w) / (t - *x))
                .sum::<Complex64>();

        let ddiff = |t: Complex64|
            zj
                .iter()
                .zip(wj.iter())
                .map(|(&x, &w)| -w / (t - x).powi(2))
                .sum::<Complex64>();

        let res: Vec<Complex64> = polzer
            .iter()
            .map(|&p| n(p) / ddiff(p))
            .collect();

        res
    }
}