use faer::{Col, Mat, Par, Spec, diag::{Diag}, dyn_stack::{MemBuffer, MemStack}, linalg::{self, qr::no_pivoting::factor::QrParams, svd::{SvdParams, bidiag::BidiagParams}}, prelude::default, traits::ComplexField, unzip, zip};
use num_complex::ComplexFloat;
use num_traits::zero;

use crate::*;

impl<T: std::fmt::Debug> Barycentric<T> {
    pub fn print_info(&self) {
        println!("Degree {:?}", &self.nodes.len());
        println!();

        println!("Error {:?}", &self.error);
        println!();

        println!("Nodes");
        for i in &self.nodes {
            println!("{:?}", i)
        };
        println!();

        println!("Values");
        for i in &self.values {
            println!("{:?}", i)
        };
        println!();

        println!("Weights");
        for i in &self.weights {
            println!("{:?}", i)
        };
        println!();

        println!("Weights * Values");
        for i in &self.weightsxvals {
            println!("{:?}", i)
        };
    }
}

pub fn calculate_weights<T>(a: &Mat<T>) -> Col<T>
where T: ComplexFloat + ComplexField, f64: From<<T as ComplexField>::Real>
{
        let mut s;
        let mut v;
        let mut w;

        let mut bdp: Spec<BidiagParams, T> = default();
        bdp.par_threshold = usize::MAX;
        bdp.config.par_threshold = usize::MAX;

        let mut qrp: Spec<QrParams, T> = default();
        qrp.blocking_threshold = usize::MAX;
        qrp.par_threshold = usize::MAX;

        let mut svdp: Spec<SvdParams, T> = default();
        svdp.recursion_threshold = usize::MAX;
        svdp.qr_ratio_threshold = f64::MAX;
        svdp.bidiag = *bdp;
        svdp.qr = *qrp;

        let (m, n) = a.shape();

        if m >= n {

            let dim: usize = Ord::min(m, n);
            s = Diag::zeros(dim);
            v = Mat::zeros(dim, dim);

            linalg::svd::svd(
                a.as_ref(),
                s.as_mut(),
                None,
                Some(v.as_mut()),
                Par::Seq,
                MemStack::new(&mut MemBuffer::new(linalg::svd::svd_scratch::<T>(
                    a.nrows(),
                    a.ncols(),
                    linalg::svd::ComputeSvdVectors::No,
                    linalg::svd::ComputeSvdVectors::Thin,
                    Par::Seq,
                    svdp,
                ))),
                svdp,
            )
            .unwrap();

            let mut s = s.into_column_vector();

            if s.iter().rev().all(|ww| *ww != zero()) {
                zip!(&mut s).for_each(|unzip!(s)| {
                    *s = s.powi(2).recip()
                });
                w = &v * &s;
                let norm: f64 = w.norm_l2().into();
                w /= norm;
            } else {
                w = v.col(v.ncols() - 1).to_owned()
            }

        } else {

            let dim: usize = Ord::min(m, n);
            let dim2: usize = Ord::max(m, n);
            s = Diag::zeros(dim);
            v = Mat::zeros(dim2, dim2);

            linalg::svd::svd(
                a.as_ref(),
                s.as_mut(),
                None,
                Some(v.as_mut()),
                Par::Seq,
                MemStack::new(&mut MemBuffer::new(linalg::svd::svd_scratch::<T>(
                    dim2,
                    dim2,
                    linalg::svd::ComputeSvdVectors::No,
                    linalg::svd::ComputeSvdVectors::Thin,
                    Par::Seq,
                    svdp,
                ))),
                svdp,
            )
            .unwrap();

            w = v.col(v.ncols() - 1).to_owned()
        }
    w
}