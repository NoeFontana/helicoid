//! The three SE(3) charts (`docs/PHASE5.md` §1.3, `0012` decision 4, `0060`).
//!
//! | Chart | `retract(δ = [φ; ρ])` | `local(Y)` |
//! |---|---|---|
//! | [`Screw`] | `X Exp(δ)` | `Log(X⁻¹ Y)` |
//! | [`Decoupled`] | `(R Exp φ, t + R ρ)` | `(Log(Rᵀ R_Y), Rᵀ (t_Y − t))` |
//! | [`WorldTranslation`] | `(R Exp φ, t + ρ)` | `(Log(Rᵀ R_Y), t_Y − t)` |
//!
//! There are two retractions, not three: `WorldTranslation` is `Decoupled` with `ρ` in the world
//! frame, `Φ(φ, ρ) = (φ, R ρ)` exactly, so it agrees with the other two to first order only if
//! `R = I` (`docs/maths/charts.md` CH.4(b), CH.5(c)). `Screw` and `Decoupled` agree to first
//! order, and their trial points differ by at most `θ/2 ‖ρ‖` (CH.7(c)). All three take a
//! [`Twist`], so the chart is a choice visible only in the type.

use crate::{Chart, Jac, LieGroup, ProductJac, RightChart, SO3Tangent, Tangent, Twist, SE3, SO3};
use helicoid_linalg::{Blend, Mat3, Matrix, Real, StridedMut};

/// The SE(3) Lie chart, `X Exp(δ)`: the group's own [`RightChart`], with `SEn3Jac` Jacobians.
pub type Screw<S> = RightChart<SE3<S>>;

/// The block-diagonal Jacobian `diag(A, B)` on a [`Twist`]: `A` on `φ`, `B` on `ρ`.
///
/// [`Decoupled`]'s and [`WorldTranslation`]'s `Jac` (CH.5), and [`SE3::chart_transition`]'s
/// result. A newtype over `ProductJac<Mat3, Mat3>` with the same blocks, dense layout and bits,
/// implementing `Jac<S, Twist<S>>` only: a second `Jac` impl on the `ProductJac` would make every
/// concrete `inverse`, `mul` or `write_dense` call ambiguous (E0283, `0060` decision 6).
#[derive(Clone, Copy, Debug)]
#[repr(transparent)]
pub struct TwistBlockJac<S>(ProductJac<Mat3<S>, Mat3<S>>);

/// The pair impl every operation but `apply` delegates to.
type Pair<S> = (SO3Tangent<S>, SO3Tangent<S>);

impl<S: Real> TwistBlockJac<S> {
    /// `diag(rotation, translation)`.
    #[inline]
    pub fn from_blocks(rotation: Mat3<S>, translation: Mat3<S>) -> Self {
        Self(ProductJac(rotation, translation))
    }
    /// The block on `φ`.
    #[inline]
    pub fn rotation_block(&self) -> Mat3<S> {
        self.0 .0
    }
    /// The block on `ρ`.
    #[inline]
    pub fn translation_block(&self) -> Mat3<S> {
        self.0 .1
    }
}

impl<S: Real> Blend<S> for TwistBlockJac<S> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self(Blend::blend(m, t.0, f.0))
    }
}

impl<S: Real> Jac<S, Twist<S>> for TwistBlockJac<S> {
    #[inline]
    fn identity() -> Self {
        Self(<ProductJac<Mat3<S>, Mat3<S>> as Jac<S, Pair<S>>>::identity())
    }
    #[inline]
    fn mul(&self, o: &Self) -> Self {
        Self(Jac::<S, Pair<S>>::mul(&self.0, &o.0))
    }
    /// Each block by `Mat3::inverse_adj`, with `Mat3`'s `Jac` domain: both blocks invertible.
    #[inline]
    fn inverse(&self) -> Self {
        Self(Jac::<S, Pair<S>>::inverse(&self.0))
    }
    #[inline]
    fn neg(&self) -> Self {
        Self(Jac::<S, Pair<S>>::neg(&self.0))
    }
    #[inline]
    fn apply(&self, t: &Twist<S>) -> Twist<S> {
        Twist {
            phi: self.0 .0 * t.phi,
            rho: [self.0 .1 * t.rho[0]],
        }
    }
    #[inline]
    fn apply_transpose(&self, t: &Twist<S>) -> Twist<S> {
        Twist {
            phi: self.0 .0.transpose() * t.phi,
            rho: [self.0 .1.transpose() * t.rho[0]],
        }
    }
    /// `Twist`'s dense order is `[φ; ρ]`, the pair's `[A; B]`: the pair's layout, unchanged.
    #[inline]
    fn write_dense(&self, out: &mut StridedMut<'_, S>) {
        Jac::<S, Pair<S>>::write_dense(&self.0, out);
    }
    #[inline]
    fn sandwich<const D: usize>(&self, cov: &Matrix<S, D, D>) -> Matrix<S, D, D> {
        const { assert!(D == <Twist<S> as Tangent<S>>::DOF) };
        Jac::<S, Pair<S>>::sandwich(&self.0, cov)
    }
}

/// `(R Exp φ, t + R ρ)`: the translation tangent in the body frame (locus-tag's `Pose::retract`).
///
/// `retract_jacobian(δ) = diag(J_r(φ), Exp(−φ))`, `local_jacobian(Y) = diag(J_r⁻¹(φ_Y), Exp(φ_Y))`,
/// `φ_Y = Log(Rᵀ R_Y)` (CH.5); `local` has `SO3::log`'s domain and `local_jacobian` `jr_inv`'s.
///
/// An iterate that linearizes in this chart, and its covariance moved to world-frame coordinates:
///
/// ```
/// use helicoid::{Chart, Decoupled, Jac, LieGroup, Manifold, Twist, WithChart, WorldTranslation};
/// use helicoid::{SE3, SO3};
/// use helicoid_linalg::{Matrix, Vector};
///
/// type Pose = WithChart<SE3<f64>, Decoupled<f64>>;
/// let x = Pose::new(SE3::from_rt(SO3::identity(), Vector([1.0, 2.0, 3.0])));
/// let chart = <Pose as Manifold<f64>>::Chart::at(&x);
/// let step = Twist { phi: Vector([0.0, 0.0, 0.1]), rho: [Vector([0.5, 0.0, 0.0])] };
/// let next: Pose = chart.retract(&step);
/// let back = chart.local(&next);
/// assert!((back.rho[0].0[0] - 0.5).abs() < 1e-15);
///
/// let cov = Matrix::<f64, 6, 6>::identity();
/// let to_world = x.0.chart_transition::<Decoupled<f64>, WorldTranslation<f64>>();
/// let _world_cov = to_world.sandwich::<6>(&cov);
/// ```
#[derive(Clone, Copy, Debug)]
pub struct Decoupled<S: Real>(SE3<S>);

/// `(R Exp φ, t + ρ)`: the translation tangent in the world frame, the chart of
/// `Product<SO3, R3>`.
///
/// `retract_jacobian(δ) = diag(J_r(φ), I)`, `local_jacobian(Y) = diag(J_r⁻¹(φ_Y), I)` (CH.5); the
/// domains are [`Decoupled`]'s.
#[derive(Clone, Copy, Debug)]
pub struct WorldTranslation<S: Real>(SE3<S>);

/// `φ` and `Rᵀ R_Y`, the rotation part both decoupled charts share.
#[inline]
fn rotation_local<S: Real>(x: &SE3<S>, y: &SE3<S>) -> (SO3Tangent<S>, SO3<S>) {
    let rel = x.rotation().inverse() * y.rotation();
    (rel.log(), rel)
}

impl<S: Real> Chart<S, SE3<S>> for Decoupled<S> {
    type Tangent = Twist<S>;
    type Jac = TwistBlockJac<S>;
    #[inline]
    fn at(base: &SE3<S>) -> Self {
        Self(*base)
    }
    #[inline]
    fn base(&self) -> SE3<S> {
        self.0
    }
    #[inline]
    fn retract(&self, delta: &Twist<S>) -> SE3<S> {
        let r = self.0.rotation();
        let phi = SO3Tangent { phi: delta.phi };
        SE3::from_rt(r.rplus(&phi), self.0.translation() + r.act(delta.rho[0]))
    }
    #[inline]
    fn local(&self, other: &SE3<S>) -> Twist<S> {
        let (phi, _) = rotation_local(&self.0, other);
        let back = self.0.rotation().inverse();
        Twist {
            phi: phi.phi,
            rho: [back.act(other.translation() - self.0.translation())],
        }
    }
    #[inline]
    fn retract_jacobian(&self, delta: &Twist<S>) -> TwistBlockJac<S> {
        let phi = SO3Tangent { phi: delta.phi };
        // `Exp(−φ)` as the transpose of `Exp(φ)`: one `Exp`, and a transpose is exact.
        let back = SO3::exp(&phi).to_matrix().transpose();
        TwistBlockJac::from_blocks(SO3::jr(&phi), back)
    }
    #[inline]
    fn local_jacobian(&self, other: &SE3<S>) -> TwistBlockJac<S> {
        // `Exp(φ_Y)` is `Rᵀ R_Y` itself, not `Exp` of its `Log`.
        let (phi, rel) = rotation_local(&self.0, other);
        TwistBlockJac::from_blocks(SO3::jr_inv(&phi), rel.to_matrix())
    }
}

impl<S: Real> Chart<S, SE3<S>> for WorldTranslation<S> {
    type Tangent = Twist<S>;
    type Jac = TwistBlockJac<S>;
    #[inline]
    fn at(base: &SE3<S>) -> Self {
        Self(*base)
    }
    #[inline]
    fn base(&self) -> SE3<S> {
        self.0
    }
    #[inline]
    fn retract(&self, delta: &Twist<S>) -> SE3<S> {
        let phi = SO3Tangent { phi: delta.phi };
        SE3::from_rt(
            self.0.rotation().rplus(&phi),
            self.0.translation() + delta.rho[0],
        )
    }
    #[inline]
    fn local(&self, other: &SE3<S>) -> Twist<S> {
        let (phi, _) = rotation_local(&self.0, other);
        Twist {
            phi: phi.phi,
            rho: [other.translation() - self.0.translation()],
        }
    }
    #[inline]
    fn retract_jacobian(&self, delta: &Twist<S>) -> TwistBlockJac<S> {
        let phi = SO3Tangent { phi: delta.phi };
        TwistBlockJac::from_blocks(SO3::jr(&phi), Matrix::identity())
    }
    #[inline]
    fn local_jacobian(&self, other: &SE3<S>) -> TwistBlockJac<S> {
        let (phi, _) = rotation_local(&self.0, other);
        TwistBlockJac::from_blocks(SO3::jr_inv(&phi), Matrix::identity())
    }
}

mod sealed {
    pub trait Sealed {}
}

impl<S: Real> sealed::Sealed for Screw<S> {}
impl<S: Real> sealed::Sealed for Decoupled<S> {}
impl<S: Real> sealed::Sealed for WorldTranslation<S> {}

/// The three SE(3) charts, and in which frame each expresses `ρ`.
///
/// **Sealed** (`0060` decision 7): `Screw`, `Decoupled` and `WorldTranslation` are the only
/// implementations, and a fourth chart is a record (`0012` decision 4). `LeftChart<SE3>` is not
/// one: its transition is `Ad_X`, not a frame change of `ρ` (CH.6).
pub trait Se3Chart<S: Real>: sealed::Sealed + Chart<S, SE3<S>, Tangent = Twist<S>> {
    /// `ρ` is a world-frame translation (`WorldTranslation`) rather than a body-frame one.
    const WORLD_TRANSLATION: bool;
}

impl<S: Real> Se3Chart<S> for Screw<S> {
    const WORLD_TRANSLATION: bool = false;
}

impl<S: Real> Se3Chart<S> for Decoupled<S> {
    const WORLD_TRANSLATION: bool = false;
}

impl<S: Real> Se3Chart<S> for WorldTranslation<S> {
    const WORLD_TRANSLATION: bool = true;
}

impl<S: Real> SE3<S> {
    /// `DΦ^{From→To}(0)` at this base (CH.4(b), CH.5(c)): `diag(I, R)` from a body-frame chart to
    /// [`WorldTranslation`], `diag(I, Rᵀ)` back, and exactly `I` otherwise, so `Screw → Decoupled`
    /// is the identity to the bit.
    ///
    /// A covariance `Σ` in `From`'s coordinates is `j.sandwich(&Σ)` in `To`'s, to first order
    /// (CH.6). `Screw` and `Decoupled` agree there but not beyond: their trial points differ by up
    /// to `θ/2 ‖ρ‖` (CH.7).
    #[inline]
    pub fn chart_transition<From: Se3Chart<S>, To: Se3Chart<S>>(&self) -> TwistBlockJac<S> {
        let translation = match (From::WORLD_TRANSLATION, To::WORLD_TRANSLATION) {
            (false, true) => self.rotation().to_matrix(),
            (true, false) => self.rotation().to_matrix().transpose(),
            _ => Matrix::identity(),
        };
        TwistBlockJac::from_blocks(Matrix::identity(), translation)
    }
}
