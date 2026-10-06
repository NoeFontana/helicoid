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
    /// The products of the dense components accumulated into `acc`, in index order.
    ///
    /// This, and not [`dot`](Tangent::dot), is the required operation (`0025`). A composite
    /// tangent — `Product<A, B>`, whose dense order is A's components then B's (`PHASE3.md` §7) —
    /// cannot build one index-order sum out of its factors' `dot`s, and stable Rust cannot size a
    /// flattening buffer from `DOF`, an associated const (`0005`); threading the accumulator gives
    /// the flat order at any nesting depth, with no buffer:
    ///
    /// ```text
    /// fn dot_acc(&self, o: &Self, acc: S) -> S {
    ///     self.1.dot_acc(&o.1, self.0.dot_acc(&o.0, acc))
    /// }
    /// ```
    fn dot_acc(&self, o: &Self, acc: S) -> S;
    /// The Euclidean inner product of the dense components, summed in index order:
    /// `dot_acc(o, +0)`.
    ///
    /// The `+0` seed is normative, not an implementation detail — `+0 + -0` is `+0`, so seeding
    /// from the first product instead would differ on a signed zero. It is what
    /// `laws::tangent_dense_order` compares against, at a recorded `f64` bound of `0`.
    fn dot(&self, o: &Self) -> S {
        self.dot_acc(o, S::zero())
    }
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
    /// `src.len() == DOF`, checked by `debug_assert!`. A release build never panics (D11) and
    /// reads a missing entry as NaN: `+0` is a valid component and would make a wrongly sized
    /// buffer indistinguishable from a zero tangent, while NaN reaches whatever the caller
    /// computes from it.
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
    /// Every impl stores entry by entry through [`StridedMut::set`], and the per-entry bounds
    /// checks are **measured to be free**: the panic branches are never taken and are laid out
    /// cold, and the loop already runs at about one store per cycle. The cost is the structural
    /// zeros, not the checks. Two designs that remove the checks were measured slower everywhere,
    /// by 53% to 326% (`0026`); do not "fix" this loop from its instruction count.
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
    /// The group geodesic `γ(x0, x1, t) = x0 · Exp(t · (x1 ⊖_R x0))` (`NUMERICS.md` §10).
    ///
    /// The provided body **is** [`reference::geodesic`](crate::reference::geodesic) and calls it,
    /// so a group that overrides this with a fast twin is measured against the one expression this
    /// returns -- `laws::geodesic_legs`'s `twin` leg, which every group carries because it costs
    /// nothing where there is no override (`PHASE4.md` §1.1, D6).
    ///
    /// `t` outside `[0, 1]` extrapolates along the same curve, and `t = 1` returns `x1` to
    /// rounding.
    ///
    /// At `t = 0` the answer is `x0` **bit for bit**, on every group that ships, with one
    /// exception: a component of `x0`'s representation that is `-0.0` can come back `+0.0`.
    /// `d.scale(0)` is `±0` per component, so the composition adds signed zeros, and a sum of
    /// zeros is negative only when every term is. The value is unchanged either way; the bit is
    /// not, and for a quaternion's `w` that sign is load-bearing, since `NUMERICS.md` §3.2 keeps
    /// `w = +0` and `Log`'s flip reads it. `sen3_tests`'s
    /// `geodesic_at_zero_is_the_left_endpoint_bit_for_bit` pins both halves.
    ///
    /// # Domain
    ///
    /// `θ(x1 ⊖_R x0) < π`. At `π` the curve is `NUMERICS.md` §3.2's function of the quaternion
    /// sign: the two preimages give geodesics `O(1)` apart, so two programs that disagree on the
    /// sign by one ulp of `θ` disagree on the answer by `O(1)`
    /// (`docs/maths/geodesics.md` GE.13(d)). Nothing panics, in debug or release.
    fn geodesic(x0: &Self, x1: &Self, t: S) -> Self {
        crate::reference::geodesic(x0, x1, t)
    }
    /// The body velocity of [`geodesic`](LieGroup::geodesic) per unit `t`, `x1 ⊖_R x0`, constant
    /// along the curve (`NUMERICS.md` §10, `docs/maths/geodesics.md` GE.2(b)). For a time step
    /// `Δt`, scale by `1/Δt`.
    fn geodesic_velocity(x0: &Self, x1: &Self) -> Self::Tangent {
        x1.rminus(x0)
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
