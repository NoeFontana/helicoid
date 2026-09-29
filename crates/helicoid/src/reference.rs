//! The reference twins: the obvious, slow implementations that the fast composites are
//! proptested against (`NUMERICS.md` §14, D6).
//!
//! A twin is the definition of *correct* for its fast counterpart, so it is public and it is not
//! optimized. Each writes its dense result into caller memory (`docs/API.md` R5: no function
//! returns a dense matrix) and keeps its scratch on the stack, sized by a const parameter `D` that
//! must equal the dimension, asserted at monomorphization (`cargo build` and `cargo test`, not
//! `cargo check`), as for [`Jac::sandwich`]:
//!
//! ```compile_fail,E0080
//! use helicoid::{reference, Jac, SEn3Jac};
//! use helicoid_linalg::StridedMut;
//! let j: SEn3Jac<f64, 1> = Jac::identity();
//! let mut out = [0.0; 49];
//! reference::sen3jac_mul::<f64, 1, 7>(&j, &j, &mut StridedMut::col_major(&mut out, 7, 7));
//! ```
//!
//! Positive control:
//!
//! ```
//! use helicoid::{reference, Jac, SEn3Jac};
//! use helicoid_linalg::StridedMut;
//! let j: SEn3Jac<f64, 1> = Jac::identity();
//! let mut out = [0.0; 36];
//! reference::sen3jac_mul::<f64, 1, 6>(&j, &j, &mut StridedMut::col_major(&mut out, 6, 6));
//! ```

use crate::dualmat::{dof, SEn3Jac};
use crate::traits::Jac;
use helicoid_linalg::{Mask, Real, StridedMut};

/// The terms added left to right from the first, `+0` when there are none.
fn sum<S: Real>(mut terms: impl Iterator<Item = S>) -> S {
    match terms.next() {
        Some(first) => terms.fold(first, |acc, t| acc + t),
        None => S::zero(),
    }
}

/// The dense image of `j` in row-major scratch, through the public [`Jac::write_dense`].
fn dense<S: Real, const N: usize, const D: usize>(j: &SEn3Jac<S, N>) -> [[S; D]; D] {
    let mut m = [[S::zero(); D]; D];
    j.write_dense(&mut StridedMut::row_major(m.as_flattened_mut(), D, D));
    m
}

/// The twin of [`SEn3Jac::mul`](Jac::mul): the dense product `dense(a) * dense(b)`, every one of
/// the `D³` products formed, structural zeros included, each entry summed left to right from
/// its first term. Writes it to the `D x D` view `out`.
///
/// # Domain
///
/// `D == 3 + 3N` (a build-time assertion) and `out` is `D x D`, checked by `debug_assert!`; a
/// smaller view panics in the strided access, the one documented panic class (D11).
pub fn sen3jac_mul<S: Real, const N: usize, const D: usize>(
    a: &SEn3Jac<S, N>,
    b: &SEn3Jac<S, N>,
    out: &mut StridedMut<'_, S>,
) {
    const { assert!(D == dof::<S, N>()) };
    debug_assert!(
        out.rows() == D && out.cols() == D,
        "sen3jac_mul: the view is not D x D"
    );
    let (x, y) = (dense::<S, N, D>(a), dense::<S, N, D>(b));
    for (r, row) in x.iter().enumerate() {
        for c in 0..D {
            out.set(r, c, sum(row.iter().zip(&y).map(|(a, yk)| *a * yk[c])));
        }
    }
}

/// The twin of [`SEn3Jac::inverse`](Jac::inverse): Gauss–Jordan elimination of `dense(j)`
/// beside the identity, which becomes the inverse. Writes it to the `D x D` view `out`.
///
/// The pivot of column `k` is the entry of largest magnitude in rows `k..`, found by comparing
/// row `k` with each later row in turn and exchanging the two where the later pivot is larger
/// (`S::select`, so no comparison reaches an `if`, `docs/API.md` R4). Its forward error is
/// `O(κ u)` for the conditioning `κ = ‖M‖_F ‖M⁻¹‖_F` of the *dense* matrix `M`, so it and the
/// structured inverse agree to a few `κ u`, not to a few `u`: the value differed by at most
/// `1.25 κ ‖M⁻¹‖_F u` in `‖·‖_F` over `10^6` seeded cases per scalar and `N = 1, 2, 3`, with
/// entries uniform on `[-1, 1)` and `κ u <= 1e-3` (`sen3jac_inverse_matches_reference`, which
/// rejects the other inputs). A nearly singular `A` is not covered by that measurement.
///
/// # Domain
///
/// `D == 3 + 3N` (a build-time assertion), `out` is `D x D` and every computed pivot is nonzero;
/// the last two are checked by `debug_assert!` (a NaN fails the pivot check). A release build
/// divides by a zero pivot and returns non-finite entries. A smaller `out` panics in the strided
/// access (D11).
///
/// ```compile_fail,E0080
/// use helicoid::{reference, Jac, SEn3Jac};
/// use helicoid_linalg::StridedMut;
/// let j: SEn3Jac<f64, 1> = Jac::identity();
/// let mut out = [0.0; 49];
/// reference::sen3jac_inverse::<f64, 1, 7>(&j, &mut StridedMut::col_major(&mut out, 7, 7));
/// ```
///
/// Positive control:
///
/// ```
/// use helicoid::{reference, Jac, SEn3Jac};
/// use helicoid_linalg::StridedMut;
/// let j: SEn3Jac<f64, 1> = Jac::identity();
/// let mut out = [0.0; 36];
/// reference::sen3jac_inverse::<f64, 1, 6>(&j, &mut StridedMut::col_major(&mut out, 6, 6));
/// ```
pub fn sen3jac_inverse<S: Real, const N: usize, const D: usize>(
    j: &SEn3Jac<S, N>,
    out: &mut StridedMut<'_, S>,
) {
    const { assert!(D == dof::<S, N>()) };
    debug_assert!(
        out.rows() == D && out.cols() == D,
        "sen3jac_inverse: the view is not D x D"
    );
    let mut m = dense::<S, N, D>(j);
    for r in 0..D {
        for c in 0..D {
            out.set(r, c, if r == c { S::one() } else { S::zero() });
        }
    }
    for k in 0..D {
        for r in k + 1..D {
            let swap = m[k][k].abs().lt(m[r][k].abs());
            let (head, tail) = m.split_at_mut(r);
            for (p, q) in head[k].iter_mut().zip(tail[0].iter_mut()) {
                (*p, *q) = (S::select(swap, *q, *p), S::select(swap, *p, *q));
            }
            for c in 0..D {
                let (p, q) = (out.get(k, c), out.get(r, c));
                out.set(k, c, S::select(swap, q, p));
                out.set(r, c, S::select(swap, p, q));
            }
        }
        let pivot = m[k][k];
        debug_assert!(
            S::zero().lt(pivot.abs()).all(),
            "sen3jac_inverse: a zero or NaN pivot"
        );
        for x in &mut m[k] {
            *x = *x / pivot;
        }
        for c in 0..D {
            out.set(k, c, out.get(k, c) / pivot);
        }
        let row = m[k];
        for (r, mr) in m.iter_mut().enumerate().filter(|&(r, _)| r != k) {
            let f = mr[k];
            for (x, y) in mr.iter_mut().zip(&row) {
                *x = *x - f * *y;
            }
            for c in 0..D {
                out.set(r, c, out.get(r, c) - f * out.get(k, c));
            }
        }
    }
}
