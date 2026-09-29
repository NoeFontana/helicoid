//! The translation group Rⁿ (`docs/PHASE3.md` §7): the trivial group that proves the traits
//! compile generically.

use crate::side::Side;
use crate::traits::{tie_dof, Jac, LieGroup, Tangent};
use core::array;
use core::ops::Mul;
use helicoid_linalg::{Blend, Mask, Matrix, Real, StridedMut, Vector};

/// An element of Rⁿ: the group law is vector addition, spelled `a * b` (`docs/API.md` R1).
///
/// `Exp` and `Log` are the identity map, so `Rn` is its own tangent space up to the wrapper
/// [`RnTangent`], and every Jacobian is `±I` or, for `ad`, `0`. There is no `Add` or `Sub`, on
/// the group or on the tangent; neither of these compiles:
///
/// ```compile_fail,E0369
/// use helicoid::Rn;
/// use helicoid_linalg::Vector;
/// let a = Rn(Vector([1.0_f64; 2]));
/// let _ = a + a;
/// ```
///
/// ```compile_fail,E0369
/// use helicoid::RnTangent;
/// use helicoid_linalg::Vector;
/// let t = RnTangent { rho: Vector([1.0_f64; 2]) };
/// let _ = t + t;
/// ```
///
/// Positive control: composition is `Mul`.
///
/// ```
/// use helicoid::Rn;
/// use helicoid_linalg::Vector;
/// let a = Rn(Vector([1.0_f64; 2]));
/// assert_eq!((a * a).0 .0, [2.0; 2]);
/// ```
#[repr(transparent)]
#[derive(Clone, Copy, Debug)]
pub struct Rn<S, const N: usize>(pub Vector<S, N>);

/// A tangent vector of [`Rn`]: `rho`, the translation, in dense order.
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct RnTangent<S, const N: usize> {
    /// The components; `write_dense` writes them in index order.
    pub rho: Vector<S, N>,
}

/// A Jacobian of [`Rn`]: the matrix `k I`, closed under `mul`, `inverse` and `neg` and holding
/// exactly the values Rⁿ needs: `I` (`Ad`, `J_r`, `J_l` and their inverses), `-I` (the second row
/// of `⊖` and `X⁻¹`) and `0` (`ad`).
#[derive(Clone, Copy, Debug)]
pub struct RnJac<S, const N: usize> {
    k: S,
}

impl<S: Real, const N: usize> Rn<S, N> {
    // The sum lives here so that clippy's `suspicious_arithmetic_impl` does not read the `+` in
    // `Mul::mul` as a typo (`Dual` does the same).
    #[inline]
    fn compose(self, o: Self) -> Self {
        Self(self.0 + o.0)
    }
}

impl<S: Real, const N: usize> Mul for Rn<S, N> {
    type Output = Self;
    #[inline]
    fn mul(self, o: Self) -> Self {
        self.compose(o)
    }
}

impl<S: Real, const N: usize> Blend<S> for Rn<S, N> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self(Vector::blend(m, t.0, f.0))
    }
}

impl<S: Real, const N: usize> Blend<S> for RnTangent<S, N> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self {
            rho: Vector::blend(m, t.rho, f.rho),
        }
    }
}

impl<S: Real, const N: usize> Blend<S> for RnJac<S, N> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self {
            k: S::select(m, t.k, f.k),
        }
    }
}

impl<S: Real, const N: usize> Tangent<S> for RnTangent<S, N> {
    const DOF: usize = N;
    #[inline]
    fn zero() -> Self {
        Self {
            rho: Vector([S::zero(); N]),
        }
    }
    #[inline]
    fn add(&self, o: &Self) -> Self {
        Self {
            rho: self.rho + o.rho,
        }
    }
    #[inline]
    fn sub(&self, o: &Self) -> Self {
        Self {
            rho: self.rho - o.rho,
        }
    }
    #[inline]
    fn neg(&self) -> Self {
        Self { rho: -self.rho }
    }
    #[inline]
    fn scale(&self, k: S) -> Self {
        Self {
            rho: self.rho.scale(k),
        }
    }
    #[inline]
    fn dot(&self, o: &Self) -> S {
        self.rho.dot(o.rho)
    }
    #[inline]
    fn write_dense(&self, out: &mut [S]) {
        debug_assert!(out.len() == N, "Tangent::write_dense: wrong length");
        for (o, v) in out.iter_mut().zip(self.rho.0) {
            *o = v;
        }
    }
    #[inline]
    fn read_dense(src: &[S]) -> Self {
        debug_assert!(src.len() == N, "Tangent::read_dense: wrong length");
        Self {
            rho: Vector(array::from_fn(|i| {
                src.get(i).copied().unwrap_or_else(S::zero)
            })),
        }
    }
}

impl<S: Real, const N: usize> RnJac<S, N> {
    // The tests build `k I` for arbitrary `k`; no public path does.
    #[inline]
    pub(crate) fn scalar(k: S) -> Self {
        Self { k }
    }
}

impl<S: Real, const N: usize> Jac<S, RnTangent<S, N>> for RnJac<S, N> {
    #[inline]
    fn identity() -> Self {
        Self::scalar(S::one())
    }
    #[inline]
    fn mul(&self, o: &Self) -> Self {
        Self::scalar(self.k * o.k)
    }
    /// # Domain
    ///
    /// `k != 0`, checked by `debug_assert!`; `ad` is `0` and has no inverse. A release build
    /// divides and returns `±inf` or NaN.
    #[inline]
    fn inverse(&self) -> Self {
        debug_assert!(S::zero().lt(self.k.abs()).all(), "RnJac::inverse: singular");
        Self::scalar(S::one() / self.k)
    }
    #[inline]
    fn neg(&self) -> Self {
        Self::scalar(-self.k)
    }
    #[inline]
    fn apply(&self, t: &RnTangent<S, N>) -> RnTangent<S, N> {
        t.scale(self.k)
    }
    // `k I` is symmetric.
    #[inline]
    fn apply_transpose(&self, t: &RnTangent<S, N>) -> RnTangent<S, N> {
        t.scale(self.k)
    }
    #[inline]
    fn write_dense(&self, out: &mut StridedMut<'_, S>) {
        debug_assert!(
            out.rows() == N && out.cols() == N,
            "Jac::write_dense: the view is not DOF x DOF"
        );
        for c in 0..N {
            for r in 0..N {
                out.set(r, c, if r == c { self.k } else { S::zero() });
            }
        }
    }
    #[inline]
    fn sandwich<const D: usize>(&self, cov: &Matrix<S, D, D>) -> Matrix<S, D, D> {
        const { assert!(D == <RnTangent<S, N> as Tangent<S>>::DOF) };
        cov.scale(self.k).scale(self.k)
    }
}

// Every row of `NUMERICS.md` §2.3 is `±I` here: `Ad = J = J⁻¹ = I` and `-J⁻¹ = -I`, the same for
// both sides, so `Sd` selects nothing in `compose_jacobians` and `inverse_jacobian`.
impl<S: Real, const N: usize> LieGroup<S> for Rn<S, N> {
    type Tangent = RnTangent<S, N>;
    type Jac = RnJac<S, N>;
    const DOF: usize = N;

    #[inline]
    fn identity() -> Self {
        const { tie_dof::<S, Self>() };
        Self(Vector([S::zero(); N]))
    }
    #[inline]
    fn inverse(&self) -> Self {
        Self(-self.0)
    }
    #[inline]
    fn exp(tau: &RnTangent<S, N>) -> Self {
        Self(tau.rho)
    }
    #[inline]
    fn log(&self) -> RnTangent<S, N> {
        RnTangent { rho: self.0 }
    }
    // `Ad = I`.
    #[inline]
    fn adjoint(&self) -> RnJac<S, N> {
        RnJac::scalar(S::one())
    }
    // The bracket of an abelian algebra vanishes.
    #[inline]
    fn ad(_tau: &RnTangent<S, N>) -> RnJac<S, N> {
        RnJac::scalar(S::zero())
    }
    #[inline]
    fn jr(_tau: &RnTangent<S, N>) -> RnJac<S, N> {
        RnJac::scalar(S::one())
    }
    #[inline]
    fn jr_inv(_tau: &RnTangent<S, N>) -> RnJac<S, N> {
        RnJac::scalar(S::one())
    }
    // `(Ad_Exp(τ)⁻¹, J_r(τ))`.
    #[inline]
    fn rplus_jacobians(&self, _tau: &RnTangent<S, N>) -> (RnJac<S, N>, RnJac<S, N>) {
        (RnJac::scalar(S::one()), RnJac::scalar(S::one()))
    }
    // `(Ad_Exp(τ), J_l(τ))`.
    #[inline]
    fn lplus_jacobians(&self, _tau: &RnTangent<S, N>) -> (RnJac<S, N>, RnJac<S, N>) {
        (RnJac::scalar(S::one()), RnJac::scalar(S::one()))
    }
    // `(J_r⁻¹, -J_l⁻¹)`.
    #[inline]
    fn rminus_jacobians(&self, _base: &Self) -> (RnJac<S, N>, RnJac<S, N>) {
        (RnJac::scalar(S::one()), RnJac::scalar(-S::one()))
    }
    // `(J_l⁻¹, -J_r⁻¹)`.
    #[inline]
    fn lminus_jacobians(&self, _base: &Self) -> (RnJac<S, N>, RnJac<S, N>) {
        (RnJac::scalar(S::one()), RnJac::scalar(-S::one()))
    }
    // Right `(Ad_Y⁻¹, I)`, left `(I, Ad_X)`.
    #[inline]
    fn compose_jacobians<Sd: Side>(&self, _rhs: &Self) -> (RnJac<S, N>, RnJac<S, N>) {
        (RnJac::scalar(S::one()), RnJac::scalar(S::one()))
    }
    // Right `-Ad_X`, left `-Ad_X⁻¹`.
    #[inline]
    fn inverse_jacobian<Sd: Side>(&self) -> RnJac<S, N> {
        RnJac::scalar(-S::one())
    }
}
