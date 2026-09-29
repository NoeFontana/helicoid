//! The three traits every group is written against (`docs/PHASE3.md` §2).

use crate::side::Side;
use core::ops::Mul;
use helicoid_linalg::{Blend, Matrix, Real, StridedMut};

/// A tangent vector: `DOF` scalars in the dense order of `NUMERICS.md` §1 (rotation first).
///
/// The operations are named methods, never `Add`/`Sub` impls (`docs/API.md` R1), and the dense
/// order is visible only through [`write_dense`](Tangent::write_dense) and
/// [`read_dense`](Tangent::read_dense) (R3).
pub trait Tangent<S: Real>: Copy + Blend<S> {
    /// The dimension. [`LieGroup::DOF`] repeats it: stable Rust cannot compute an array length
    /// from an associated const, so the two are tied by a `const` assertion (`0005`).
    const DOF: usize;
    /// The zero tangent.
    fn zero() -> Self;
    /// `self + o`.
    fn add(&self, o: &Self) -> Self;
    /// `self - o`.
    fn sub(&self, o: &Self) -> Self;
    /// `-self`.
    fn neg(&self) -> Self;
    /// `k * self`.
    fn scale(&self, k: S) -> Self;
    /// The Euclidean inner product of the dense components, summed in index order.
    fn dot(&self, o: &Self) -> S;
    /// Writes the components in the order of `NUMERICS.md` §1.
    ///
    /// # Domain
    ///
    /// `out.len() == DOF`, checked by `debug_assert!`. A release build never panics and writes
    /// the first `min(out.len(), DOF)` entries.
    fn write_dense(&self, out: &mut [S]);
    /// The tangent whose dense components are `src`.
    ///
    /// # Domain
    ///
    /// `src.len() == DOF`, checked by `debug_assert!`. A release build never panics and reads a
    /// missing entry as `+0`.
    fn read_dense(src: &[S]) -> Self;
}

/// A Jacobian of a group operation, expressed on tangents of type `T` (`docs/API.md` R5).
///
/// A structured type closed under [`mul`](Jac::mul), [`inverse`](Jac::inverse) and
/// [`neg`](Jac::neg) (`0005`). A dense matrix leaves it only through
/// [`write_dense`](Jac::write_dense) and [`sandwich`](Jac::sandwich).
pub trait Jac<S: Real, T: Tangent<S>>: Copy + Blend<S> {
    /// The identity matrix.
    fn identity() -> Self;
    /// The matrix product `self * o`.
    fn mul(&self, o: &Self) -> Self;
    /// The exact algebraic inverse, not a numerical one (`0005`).
    ///
    /// # Domain
    ///
    /// `self` is invertible. Each impl states what it checks; a release build never panics.
    fn inverse(&self) -> Self;
    /// `-self`.
    fn neg(&self) -> Self;
    /// `self * t`.
    fn apply(&self, t: &T) -> T;
    /// `selfᵀ * t`.
    fn apply_transpose(&self, t: &T) -> T;
    /// Writes the dense `DOF x DOF` matrix into `out`: row `r` and column `c` are indices of the
    /// dense order of `T`, and every entry is written, structural zeros included.
    ///
    /// # Domain
    ///
    /// `out` is `DOF x DOF`, checked by `debug_assert!`. A smaller view panics in the strided
    /// access, the one documented panic class (D11).
    fn write_dense(&self, out: &mut StridedMut<'_, S>);
    /// `self * cov * selfᵀ` for a dense covariance.
    ///
    /// Every impl starts with `const { assert!(D == T::DOF) }`, so a covariance of the wrong size
    /// fails to build, at monomorphization (`cargo build` and `cargo test`, not `cargo check`):
    ///
    /// ```compile_fail,E0080
    /// use helicoid::{Jac, RnJac};
    /// use helicoid_linalg::Matrix;
    /// let j: RnJac<f64, 3> = Jac::identity();
    /// let _ = j.sandwich::<4>(&Matrix::<f64, 4, 4>::identity());
    /// ```
    ///
    /// Positive control, so the failure above comes from the assertion and not from an import:
    ///
    /// ```
    /// use helicoid::{Jac, RnJac};
    /// use helicoid_linalg::Matrix;
    /// let j: RnJac<f64, 3> = Jac::identity();
    /// let _ = j.sandwich::<3>(&Matrix::<f64, 3, 3>::identity());
    /// ```
    fn sandwich<const D: usize>(&self, cov: &Matrix<S, D, D>) -> Matrix<S, D, D>;
}

/// A Lie group over the scalar `S`: composition is `Mul`, tangents are rotation-first
/// (`docs/API.md` R1, `0002`).
///
/// Perturbations are always spelled with their side: `rplus`, `lplus`, `rminus`, `lminus`
/// (`X ⊕_R τ = X · Exp(τ)`, `Y ⊖_R X = Log(X⁻¹ · Y)`, and the left forms). Each `*_jacobians`
/// method returns `(∂/∂first, ∂/∂second)` in its side's convention, the rows of `NUMERICS.md`
/// §2.3. A provided method is overridden only where a closed form is cheaper, and the provided
/// body then becomes the reference twin (`NUMERICS.md` §14).
pub trait LieGroup<S: Real>: Copy + Blend<S> + Mul<Output = Self> {
    /// The tangent space, rotation-first.
    type Tangent: Tangent<S>;
    /// The structured Jacobian type (`0005`).
    type Jac: Jac<S, Self::Tangent>;
    /// The dimension, equal to `Self::Tangent::DOF`. Every impl asserts the tie at the top of
    /// [`identity`](LieGroup::identity), `const { assert!(Self::DOF == <Self::Tangent as
    /// Tangent<S>>::DOF) }`, which fires when `identity` is instantiated: a mismatch fails
    /// `cargo build` and `cargo test`, not `cargo check`. The test-side laws assert it as well.
    const DOF: usize;
    /// The identity element.
    fn identity() -> Self;
    /// The group inverse.
    fn inverse(&self) -> Self;
    /// `Exp(τ)`.
    fn exp(tau: &Self::Tangent) -> Self;
    /// `Log(self)`, the inverse of `Exp` on its canonical branch.
    fn log(&self) -> Self::Tangent;
    /// `Ad_self`: `self · Exp(τ) · self⁻¹ = Exp(Ad_self τ)`.
    fn adjoint(&self) -> Self::Jac;
    /// `ad_τ`: `ad_τ σ = [τ^, σ^]^vee`.
    fn ad(tau: &Self::Tangent) -> Self::Jac;
    /// The right Jacobian of `Exp`: `Exp(τ + δ) ≈ Exp(τ) · Exp(J_r δ)`.
    fn jr(tau: &Self::Tangent) -> Self::Jac;
    /// `J_r⁻¹`.
    ///
    /// # Domain
    ///
    /// Each impl states its domain (`NUMERICS.md` §12).
    fn jr_inv(tau: &Self::Tangent) -> Self::Jac;

    /// The left Jacobian of `Exp`: `J_l(τ) = J_r(-τ)`.
    fn jl(tau: &Self::Tangent) -> Self::Jac {
        Self::jr(&tau.neg())
    }
    /// `J_l⁻¹`; see [`jr_inv`](LieGroup::jr_inv) for the domain.
    fn jl_inv(tau: &Self::Tangent) -> Self::Jac {
        Self::jr_inv(&tau.neg())
    }
    /// `self ⊕_R τ = self · Exp(τ)`.
    fn rplus(&self, tau: &Self::Tangent) -> Self {
        *self * Self::exp(tau)
    }
    /// `self ⊕_L τ = Exp(τ) · self`.
    fn lplus(&self, tau: &Self::Tangent) -> Self {
        Self::exp(tau) * *self
    }
    /// `self ⊖_R base = Log(base⁻¹ · self)`.
    fn rminus(&self, base: &Self) -> Self::Tangent {
        (base.inverse() * *self).log()
    }
    /// `self ⊖_L base = Log(self · base⁻¹)`.
    fn lminus(&self, base: &Self) -> Self::Tangent {
        (*self * base.inverse()).log()
    }
    /// `(∂(self ⊕_R τ)/∂self, ∂(self ⊕_R τ)/∂τ)`, right convention.
    fn rplus_jacobians(&self, tau: &Self::Tangent) -> (Self::Jac, Self::Jac);
    /// `(∂(self ⊕_L τ)/∂self, ∂(self ⊕_L τ)/∂τ)`, left convention.
    fn lplus_jacobians(&self, tau: &Self::Tangent) -> (Self::Jac, Self::Jac);
    /// `(∂(self ⊖_R base)/∂self, ∂(self ⊖_R base)/∂base)`, right convention.
    fn rminus_jacobians(&self, base: &Self) -> (Self::Jac, Self::Jac);
    /// `(∂(self ⊖_L base)/∂self, ∂(self ⊖_L base)/∂base)`, left convention.
    fn lminus_jacobians(&self, base: &Self) -> (Self::Jac, Self::Jac);
    /// `(∂(self · rhs)/∂self, ∂(self · rhs)/∂rhs)` in the convention of `Sd`.
    fn compose_jacobians<Sd: Side>(&self, rhs: &Self) -> (Self::Jac, Self::Jac);
    /// `∂(self⁻¹)/∂self` in the convention of `Sd`.
    fn inverse_jacobian<Sd: Side>(&self) -> Self::Jac;
}

/// Ties `G::DOF` to `G::Tangent::DOF`; `const { tie_dof::<S, Self>() }` at the top of an impl's
/// `identity` fails the build when they differ.
pub(crate) const fn tie_dof<S: Real, G: LieGroup<S>>() {
    assert!(
        G::DOF == <G::Tangent as Tangent<S>>::DOF,
        "LieGroup::DOF must equal Tangent::DOF"
    );
}
