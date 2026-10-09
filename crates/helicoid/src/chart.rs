//! Charts: a retraction and its local inverse, frozen at a base (`docs/PHASE5.md` §1, `0012`,
//! `0060`).
//!
//! A chart at `X` is the pair `retract: δ ↦ Y` and `local: Y ↦ δ`, with `local(retract(δ)) = δ`
//! for `δ` in `U₀ = {θ < π}` (`docs/maths/charts.md` CH.1). Its two Jacobians read the perturbed
//! point **in the chart frozen at that point**: `retract_jacobian(δ)` is
//! `∂ local_Y(retract_X(δ + η)) / ∂η` at `η = 0`, `Y = retract_X(δ)`, and `local_jacobian(Y)` is
//! `∂ local_X(retract_Y(η)) / ∂η`, so `local_jacobian(retract(δ)) = retract_jacobian(δ)⁻¹`
//! (CH.2(b)).
//!
//! Every `LieGroup` is a [`Manifold`] whose chart is [`RightChart`]. There is one impl per group
//! and no blanket over `LieGroup`: beside [`WithChart`]'s impl a blanket is E0119, since a
//! downstream crate may implement `LieGroup<ItsReal>` for a `WithChart` (`0060` decision 4).

use crate::{LieGroup, Product, Rn, SEn3, Tangent, SO3};
use core::{fmt, marker::PhantomData};
use helicoid_linalg::{Blend, Real};

/// A chart of `M`, frozen at a base by [`at`](Chart::at) (`docs/PHASE5.md` §1.1).
///
/// `retract`, `local` and their Jacobians use only what `at` froze: a chart never recomputes its
/// frame from the point it is applied to (`0012` decision 1).
pub trait Chart<S: Real, M>: Copy {
    /// The chart's coordinates.
    type Tangent: Tangent<S>;
    /// The Jacobians of [`retract`](Chart::retract) and [`local`](Chart::local), on those
    /// coordinates.
    type Jac: crate::Jac<S, Self::Tangent>;
    /// The chart at `base`.
    fn at(base: &M) -> Self;
    /// The base it was frozen at.
    fn base(&self) -> M;
    /// The point at coordinates `delta`.
    fn retract(&self, delta: &Self::Tangent) -> M;
    /// The coordinates of `other`.
    ///
    /// # Domain
    ///
    /// Each chart states its own; `local(retract(δ)) = δ` holds on `θ < π` (CH.2).
    fn local(&self, other: &M) -> Self::Tangent;
    /// `∂ local_Y(retract(δ + η)) / ∂η` at `η = 0`, `Y = retract(δ)` (CH.1). A formula defined for
    /// every `δ`, invertible on `θ < π`.
    fn retract_jacobian(&self, delta: &Self::Tangent) -> Self::Jac;
    /// `∂ local(retract_Y(η)) / ∂η` at `η = 0`, `Y = other` (CH.1): the inverse of
    /// `retract_jacobian(local(other))` (CH.2(b)).
    ///
    /// # Domain
    ///
    /// That of [`local`](Chart::local), and each chart states what its inverse Jacobian needs.
    fn local_jacobian(&self, other: &M) -> Self::Jac;
}

/// A space a solver linearizes in its default chart (`docs/PHASE5.md` §1.1).
///
/// With [`LieGroup`] also in scope, a bare `SE3::<f64>::DOF` is ambiguous (E0034): write
/// `<SE3<f64> as LieGroup<f64>>::DOF`. The two are equal on every group.
pub trait Manifold<S: Real>: Copy + Blend<S> {
    /// The coordinates of [`Chart`](Manifold::Chart).
    type Tangent: Tangent<S>;
    /// The chart a solver linearizes in by default.
    type Chart: Chart<S, Self, Tangent = Self::Tangent>;
    /// The dimension, `Tangent::DOF`.
    const DOF: usize;
}

/// The right chart of a group: `retract(δ) = X Exp(δ)`, `local(Y) = Y ⊖_R X` (CH.3).
///
/// `retract_jacobian(δ) = J_r(δ)`; `local_jacobian(Y) = J_r⁻¹(Y ⊖_R X)`, with `jr_inv`'s domain.
#[derive(Clone, Copy, Debug)]
pub struct RightChart<G>(G);

/// The left chart of a group: `retract(δ) = Exp(δ) X`, `local(Y) = Y ⊖_L X` (CH.3).
///
/// `retract_jacobian(δ) = J_l(δ)`; `local_jacobian(Y) = J_l⁻¹(Y ⊖_L X)`, with `jl_inv`'s domain.
/// It is the right chart after the linear map `Ad_X` (CH.6).
#[derive(Clone, Copy, Debug)]
pub struct LeftChart<G>(G);

impl<S: Real, G: LieGroup<S>> Chart<S, G> for RightChart<G> {
    type Tangent = G::Tangent;
    type Jac = G::Jac;
    #[inline]
    fn at(base: &G) -> Self {
        Self(*base)
    }
    #[inline]
    fn base(&self) -> G {
        self.0
    }
    #[inline]
    fn retract(&self, delta: &G::Tangent) -> G {
        self.0.rplus(delta)
    }
    #[inline]
    fn local(&self, other: &G) -> G::Tangent {
        other.rminus(&self.0)
    }
    #[inline]
    fn retract_jacobian(&self, delta: &G::Tangent) -> G::Jac {
        G::jr(delta)
    }
    #[inline]
    fn local_jacobian(&self, other: &G) -> G::Jac {
        G::jr_inv(&self.local(other))
    }
}

impl<S: Real, G: LieGroup<S>> Chart<S, G> for LeftChart<G> {
    type Tangent = G::Tangent;
    type Jac = G::Jac;
    #[inline]
    fn at(base: &G) -> Self {
        Self(*base)
    }
    #[inline]
    fn base(&self) -> G {
        self.0
    }
    #[inline]
    fn retract(&self, delta: &G::Tangent) -> G {
        self.0.lplus(delta)
    }
    #[inline]
    fn local(&self, other: &G) -> G::Tangent {
        other.lminus(&self.0)
    }
    #[inline]
    fn retract_jacobian(&self, delta: &G::Tangent) -> G::Jac {
        G::jl(delta)
    }
    #[inline]
    fn local_jacobian(&self, other: &G) -> G::Jac {
        G::jl_inv(&self.local(other))
    }
}

/// `M` linearized in the chart `C` instead of its default, so the choice is in the type
/// (`0012` decision 3): locus-tag's iterate is `WithChart<SE3<f64>, Decoupled<f64>>`.
#[repr(transparent)]
pub struct WithChart<M, C>(pub M, PhantomData<C>);

impl<M, C> WithChart<M, C> {
    /// `m`, linearized in `C`.
    #[inline]
    pub const fn new(m: M) -> Self {
        Self(m, PhantomData)
    }
}

// Written on `M` alone: a derive would also ask `C` for each trait, and a chart need not be `Debug`.
impl<M: Clone, C> Clone for WithChart<M, C> {
    #[inline]
    fn clone(&self) -> Self {
        Self::new(self.0.clone())
    }
}

impl<M: Copy, C> Copy for WithChart<M, C> {}

impl<M: fmt::Debug, C> fmt::Debug for WithChart<M, C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_tuple("WithChart").field(&self.0).finish()
    }
}

impl<S: Real, M: Blend<S>, C> Blend<S> for WithChart<M, C> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self::new(M::blend(m, t.0, f.0))
    }
}

/// [`WithChart`]'s chart: `C`, wrapping and unwrapping the point.
///
/// A second `Chart` impl on `C` itself would give every concrete `retract` call two candidates
/// (`0060` decision 5); this is the only `Chart` over a `WithChart`.
#[derive(Clone, Copy, Debug)]
#[repr(transparent)]
pub struct Lifted<C>(C);

impl<S: Real, M: Copy, C: Chart<S, M>> Chart<S, WithChart<M, C>> for Lifted<C> {
    type Tangent = C::Tangent;
    type Jac = C::Jac;
    #[inline]
    fn at(base: &WithChart<M, C>) -> Self {
        Self(C::at(&base.0))
    }
    #[inline]
    fn base(&self) -> WithChart<M, C> {
        WithChart::new(self.0.base())
    }
    #[inline]
    fn retract(&self, delta: &C::Tangent) -> WithChart<M, C> {
        WithChart::new(self.0.retract(delta))
    }
    #[inline]
    fn local(&self, other: &WithChart<M, C>) -> C::Tangent {
        self.0.local(&other.0)
    }
    #[inline]
    fn retract_jacobian(&self, delta: &C::Tangent) -> C::Jac {
        self.0.retract_jacobian(delta)
    }
    #[inline]
    fn local_jacobian(&self, other: &WithChart<M, C>) -> C::Jac {
        self.0.local_jacobian(&other.0)
    }
}

impl<S: Real, M: Copy + Blend<S>, C: Chart<S, M>> Manifold<S> for WithChart<M, C> {
    type Tangent = C::Tangent;
    type Chart = Lifted<C>;
    const DOF: usize = <C::Tangent as Tangent<S>>::DOF;
}

// One `Manifold` impl per group, `Chart = RightChart<Self>`, written out (`0060` decision 4): see
// the module header for why there is no blanket. A group added later adds its impl here.

impl<S: Real> Manifold<S> for SO3<S> {
    type Tangent = <Self as LieGroup<S>>::Tangent;
    type Chart = RightChart<Self>;
    const DOF: usize = <Self as LieGroup<S>>::DOF;
}

impl<S: Real, const N: usize> Manifold<S> for SEn3<S, N> {
    type Tangent = <Self as LieGroup<S>>::Tangent;
    type Chart = RightChart<Self>;
    const DOF: usize = <Self as LieGroup<S>>::DOF;
}

impl<S: Real, const N: usize> Manifold<S> for Rn<S, N> {
    type Tangent = <Self as LieGroup<S>>::Tangent;
    type Chart = RightChart<Self>;
    const DOF: usize = <Self as LieGroup<S>>::DOF;
}

impl<S: Real, A: LieGroup<S>, B: LieGroup<S>> Manifold<S> for Product<A, B> {
    type Tangent = <Self as LieGroup<S>>::Tangent;
    type Chart = RightChart<Self>;
    const DOF: usize = <Self as LieGroup<S>>::DOF;
}
