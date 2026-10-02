//! A non-abelian group under the generic laws of [`crate::laws`]: the Heisenberg group in
//! exponential coordinates, test-only.
//!
//! `Rn` cannot tell the sides apart (`Ad = J = I`) and has only symmetric Jacobians, so a swapped
//! `Left`/`Right` wiring, a wrong `lplus`, `lminus`, `jl` or `jl_inv` body and a transposed `Jac`
//! pass its tests. Here `Exp` and `Log` are the identity on coordinates,
//! `X Y = X + Y + [X, Y] / 2` with `[X, Y] = (0, 0, x1 y2 - x2 y1)` and `ad² = 0`, so
//! `Ad = I + ad`, `J_r = I - ad / 2`, `J_l = I + ad / 2` and every row of `NUMERICS.md` §2.3 is
//! a polynomial. The closed forms below were checked against the perturbation definitions of those
//! rows on the 3x3 unipotent matrices (matrix `exp` and `log` series, sympy). The group is 2-step
//! nilpotent: it checks no series term past the first, and its `Log` has no branch.
//!
//! The bounds are recorded as described at [`Bounds`]. The worst errors, in `u`: `f64` and
//! `Dual<f64, 3>` 4.00 for `group_axioms`, 3.85 for `adjoint_identity`, 3.97 for `plus_minus`, 5.89
//! for `jac_dense_order`, 5.55 for `sandwich_matches_dense` and 0 for the rest; `f32` 4.08, 3.53,
//! 3.95, 5.48, 4.09, 1.61 for `tangent_dense_order` and 0 for the rest.

use crate::laws::{dense_bits, laws_for, Bounds, Sample};
use crate::{Jac, LieGroup, Right, RnTangent, Side, Tangent};
use core::any::TypeId;
use core::array;
use core::ops::Mul;
use helicoid_linalg::{Blend, Mat3, Matrix, Real, StridedMut, Vector};
use std::vec::Vec;

/// `(x1, x2, x3)`: the exponential coordinates of `[[1, x1, x3 + x1 x2 / 2], [0, 1, x2], [0, 0, 1]]`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Heis<S>(pub(crate) Vector<S, 3>);

/// A dense `3 x 3` Jacobian, inverted by its adjugate; it states no domain, since `probe` inverts
/// NaN matrices.
#[derive(Clone, Copy, Debug)]
pub(crate) struct HJac<S>(pub(crate) Mat3<S>);

impl<S: Real> Heis<S> {
    // The cross term lives here so that clippy's `suspicious_arithmetic_impl` does not read the `-`
    // in `Mul::mul` as a typo (as in `Rn`).
    fn compose(self, o: Self) -> Self {
        let (x, y) = (self.0 .0, o.0 .0);
        let cross = S::lit(0.5) * (x[0] * y[1] - x[1] * y[0]);
        Self(Vector([x[0] + y[0], x[1] + y[1], x[2] + y[2] + cross]))
    }
}

impl<S: Real> Mul for Heis<S> {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        self.compose(o)
    }
}

impl<S: Real> Blend<S> for Heis<S> {
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self(Vector::blend(m, t.0, f.0))
    }
}

impl<S: Real> Blend<S> for HJac<S> {
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self(Matrix::blend(m, t.0, f.0))
    }
}

/// `ad_t`: `ad_t s = (0, 0, t1 s2 - t2 s1)`.
fn ad_matrix<S: Real>(t: &Vector<S, 3>) -> Mat3<S> {
    let z = Vector([S::zero(); 3]);
    Matrix::from_rows([z, z, Vector([-t.0[1], t.0[0], S::zero()])])
}

/// `I + k ad_t`.
fn jac<S: Real>(t: &Vector<S, 3>, k: S) -> HJac<S> {
    HJac(Matrix::identity() + ad_matrix(t).scale(k))
}

impl<S: Real> Jac<S, RnTangent<S, 3>> for HJac<S> {
    fn identity() -> Self {
        Self(Matrix::identity())
    }
    fn mul(&self, o: &Self) -> Self {
        Self(self.0 * o.0)
    }
    fn inverse(&self) -> Self {
        Self(self.0.inverse_adj().0)
    }
    fn neg(&self) -> Self {
        Self(-self.0)
    }
    fn apply(&self, t: &RnTangent<S, 3>) -> RnTangent<S, 3> {
        RnTangent {
            rho: self.0 * t.rho,
        }
    }
    fn apply_transpose(&self, t: &RnTangent<S, 3>) -> RnTangent<S, 3> {
        RnTangent {
            rho: self.0.transpose() * t.rho,
        }
    }
    fn write_dense(&self, out: &mut StridedMut<'_, S>) {
        debug_assert!(
            out.rows() == 3 && out.cols() == 3,
            "Jac::write_dense: the view is not DOF x DOF"
        );
        for c in 0..3 {
            for r in 0..3 {
                out.set(r, c, self.0.get(r, c));
            }
        }
    }
    fn sandwich<const D: usize>(&self, cov: &Matrix<S, D, D>) -> Matrix<S, D, D> {
        const { assert!(D == <RnTangent<S, 3> as Tangent<S>>::DOF) };
        let m = Mat3::from_cols(array::from_fn(|c| {
            Vector(array::from_fn(|r| cov.get(r, c)))
        }));
        let s = self.0 * m * self.0.transpose();
        Matrix::from_cols(array::from_fn(|c| Vector(array::from_fn(|r| s.get(r, c)))))
    }
}

/// `Sd` is `Right`: `Side` has no method that selects a side row, so this group compares types.
fn is_right<Sd: Side>() -> bool {
    TypeId::of::<Sd>() == TypeId::of::<Right>()
}

// The provided `jl`, `jl_inv`, `rplus`, `lplus`, `rminus` and `lminus` are left alone, so the laws
// check their bodies; the rows are the closed forms of `NUMERICS.md` §2.3.
impl<S: Real> LieGroup<S> for Heis<S> {
    type Tangent = RnTangent<S, 3>;
    type Jac = HJac<S>;
    const DOF: usize = 3;

    fn identity() -> Self {
        const { crate::traits::tie_dof::<S, Self>() };
        Self(Vector([S::zero(); 3]))
    }
    fn inverse(&self) -> Self {
        Self(-self.0)
    }
    fn exp(tau: &RnTangent<S, 3>) -> Self {
        Self(tau.rho)
    }
    fn log(&self) -> RnTangent<S, 3> {
        RnTangent { rho: self.0 }
    }
    fn adjoint(&self) -> HJac<S> {
        jac(&self.0, S::one())
    }
    fn ad(tau: &RnTangent<S, 3>) -> HJac<S> {
        HJac(ad_matrix(&tau.rho))
    }
    fn jr(tau: &RnTangent<S, 3>) -> HJac<S> {
        jac(&tau.rho, -S::lit(0.5))
    }
    fn jr_inv(tau: &RnTangent<S, 3>) -> HJac<S> {
        jac(&tau.rho, S::lit(0.5))
    }
    fn rplus_jacobians(&self, tau: &RnTangent<S, 3>) -> (HJac<S>, HJac<S>) {
        (jac(&tau.rho, -S::one()), jac(&tau.rho, -S::lit(0.5)))
    }
    fn lplus_jacobians(&self, tau: &RnTangent<S, 3>) -> (HJac<S>, HJac<S>) {
        (jac(&tau.rho, S::one()), jac(&tau.rho, S::lit(0.5)))
    }
    fn rminus_jacobians(&self, base: &Self) -> (HJac<S>, HJac<S>) {
        let tau = (base.inverse() * *self).0;
        (jac(&tau, S::lit(0.5)), jac(&tau, -S::lit(0.5)).neg())
    }
    fn lminus_jacobians(&self, base: &Self) -> (HJac<S>, HJac<S>) {
        let tau = (*self * base.inverse()).0;
        (jac(&tau, -S::lit(0.5)), jac(&tau, S::lit(0.5)).neg())
    }
    fn compose_jacobians<Sd: Side>(&self, rhs: &Self) -> (HJac<S>, HJac<S>) {
        if is_right::<Sd>() {
            (jac(&rhs.0, -S::one()), HJac::identity())
        } else {
            (HJac::identity(), jac(&self.0, S::one()))
        }
    }
    fn inverse_jacobian<Sd: Side>(&self) -> HJac<S> {
        let k = if is_right::<Sd>() {
            S::one()
        } else {
            -S::one()
        };
        jac(&self.0, k).neg()
    }
}

/// `I + D M / 4` with `M[r][c] = v[(c - r) mod 3]` the circulant of the sample and
/// `D = diag(1, 1/2, 1/4)`: full, non-symmetric and strictly diagonally dominant, so invertible
/// with a condition number below 7.
///
/// The row scaling is what makes the family **non-commuting**, and that is the whole point of this
/// group: `jac_dense_order` pins `mul`'s operand order by comparing `a.mul(b)` against the dense
/// product in that order, which detects nothing when `a b = b a`. A circulant commutes with every
/// circulant, and `I` plus a circulant is one, so with `M` alone a `mul` that multiplies its
/// operands the wrong way round passed every law; `D M D' M'` and `D M' D M` differ.
///
/// `lane` offsets the `Dual` variable each sample becomes. A fixture that holds this block beside
/// others gives each one the lanes of its own place in the dense order, so no two blocks claim a
/// lane while holding different values.
pub(crate) fn heis_jac_at<S: Sample>(v: &[f64; 3], lane: usize) -> HJac<S> {
    let row = |r: usize| {
        let k = S::lit(0.25 / f64::from(1_u8 << r));
        Vector(array::from_fn(|c| {
            let i = (c + 3 - r) % 3;
            k * S::sample(v[i], lane + i)
        }))
    };
    HJac(Matrix::identity() + Matrix::from_rows([row(0), row(1), row(2)]))
}

/// [`heis_jac_at`] at lane `0`, the one-argument form `laws_for!` calls.
fn heis_jac<S: Sample>(v: &[f64; 3]) -> HJac<S> {
    heis_jac_at(v, 0)
}

const F64: Bounds = Bounds {
    axioms: 8.0,
    exp_log: 0.0,
    adjoint: 8.0,
    jl_ad_jr: 0.0,
    plus_minus: 8.0,
    rows: 0.0,
    ad: 0.0,
    sides: 0.0,
    tangent_order: 0.0,
    jac_order: 12.0,
    sandwich: 12.0,
};
const F32: Bounds = Bounds {
    axioms: 9.0,
    tangent_order: 4.0,
    jac_order: 11.0,
    sandwich: 9.0,
    ..F64
};

laws_for!(heis, Heis, heis_jac, 3, F64, F32);

fn bits(x: &[f64]) -> Vec<u64> {
    x.iter().map(|v| v.to_bits()).collect()
}

// No law calls `blend`, and neither type has a `branch`/`select` that would reach it, so without
// this the two impls are unexecuted and swapping their arms changes nothing.
#[test]
fn blend_selects_each_value_type() {
    let (a, b) = (Heis(Vector([1.0; 3])), Heis(Vector([2.0; 3])));
    let (ja, jb) = (heis_jac::<f64>(&[0.5; 3]), heis_jac::<f64>(&[-0.5; 3]));
    for (m, x, j) in [(true, a, ja), (false, b, jb)] {
        assert_eq!(bits(&Heis::blend(m, a, b).0 .0), bits(&x.0 .0));
        assert_eq!(
            dense_bits::<_, _, 3>(&HJac::blend(m, ja, jb)),
            dense_bits::<_, _, 3>(&j)
        );
    }
}

fn cols(m: [[f64; 3]; 3]) -> [[u64; 3]; 3] {
    m.map(|c| c.map(f64::to_bits))
}

// Hand values, so that the group is what its header says and the sides really differ.
#[test]
fn the_sides_differ() {
    let (x, y) = (Heis(Vector([1.0, 2.0, 3.0])), Heis(Vector([4.0, 5.0, 6.0])));
    let tau = y.log();
    assert_eq!(bits(&(x * y).0 .0), bits(&[5.0, 7.0, 7.5]));
    assert_eq!(bits(&x.rplus(&tau).0 .0), bits(&[5.0, 7.0, 7.5]));
    assert_eq!(bits(&x.lplus(&tau).0 .0), bits(&[5.0, 7.0, 10.5]));
    assert_eq!(bits(&y.rminus(&x).rho.0), bits(&[3.0, 3.0, 4.5]));
    assert_eq!(bits(&y.lminus(&x).rho.0), bits(&[3.0, 3.0, 1.5]));
    let sigma = RnTangent {
        rho: Vector([1.0, 2.0, 3.0]),
    };
    assert_eq!(
        bits(&Heis::<f64>::ad(&tau).apply(&sigma).rho.0),
        bits(&[0.0, 0.0, 3.0])
    );
    // Columns of `ad_tau`, `J_l(tau)` and `J_r(tau)`: not symmetric.
    let ad = cols([[0.0, 0.0, -5.0], [0.0, 0.0, 4.0], [0.0; 3]]);
    let jl = cols([[1.0, 0.0, -2.5], [0.0, 1.0, 2.0], [0.0, 0.0, 1.0]]);
    let jr = cols([[1.0, 0.0, 2.5], [0.0, 1.0, -2.0], [0.0, 0.0, 1.0]]);
    assert_eq!(dense_bits::<_, _, 3>(&Heis::<f64>::ad(&tau)), ad);
    assert_eq!(dense_bits::<_, _, 3>(&Heis::<f64>::jl(&tau)), jl);
    assert_eq!(dense_bits::<_, _, 3>(&Heis::<f64>::jr(&tau)), jr);
}
