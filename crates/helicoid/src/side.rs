//! Which end a perturbation multiplies (`docs/PHASE3.md` §2, `docs/API.md` R1).

use crate::traits::LieGroup;
use helicoid_linalg::Real;

/// The end of a product a perturbation multiplies: [`Right`] (`X · Exp(τ)`, the default) or
/// [`Left`] (`Exp(τ) · X`). Generic code takes `Sd: Side`, never a flag, and each side's Jacobians
/// are expressed in that side's convention (`NUMERICS.md` §2.3).
pub trait Side: Copy + 'static {
    /// `x ⊕ τ`: [`LieGroup::rplus`] for `Right`, [`LieGroup::lplus`] for `Left`.
    fn plus<S: Real, G: LieGroup<S>>(x: &G, tau: &G::Tangent) -> G;
    /// `y ⊖ x`: [`LieGroup::rminus`] for `Right`, [`LieGroup::lminus`] for `Left`.
    fn minus<S: Real, G: LieGroup<S>>(y: &G, x: &G) -> G::Tangent;
    /// `(∂(x ⊕ τ)/∂x, ∂(x ⊕ τ)/∂τ)`: `rplus_jacobians` or `lplus_jacobians`.
    fn plus_jacobians<S: Real, G: LieGroup<S>>(x: &G, tau: &G::Tangent) -> (G::Jac, G::Jac);
    /// `(∂(y ⊖ x)/∂y, ∂(y ⊖ x)/∂x)`: `rminus_jacobians` or `lminus_jacobians`.
    fn minus_jacobians<S: Real, G: LieGroup<S>>(y: &G, x: &G) -> (G::Jac, G::Jac);
}

/// Right perturbation, `X · Exp(τ)`: the default side (`0002`).
#[derive(Clone, Copy, Debug)]
pub struct Right;

/// Left perturbation, `Exp(τ) · X`.
#[derive(Clone, Copy, Debug)]
pub struct Left;

impl Side for Right {
    fn plus<S: Real, G: LieGroup<S>>(x: &G, tau: &G::Tangent) -> G {
        x.rplus(tau)
    }
    fn minus<S: Real, G: LieGroup<S>>(y: &G, x: &G) -> G::Tangent {
        y.rminus(x)
    }
    fn plus_jacobians<S: Real, G: LieGroup<S>>(x: &G, tau: &G::Tangent) -> (G::Jac, G::Jac) {
        x.rplus_jacobians(tau)
    }
    fn minus_jacobians<S: Real, G: LieGroup<S>>(y: &G, x: &G) -> (G::Jac, G::Jac) {
        y.rminus_jacobians(x)
    }
}

impl Side for Left {
    fn plus<S: Real, G: LieGroup<S>>(x: &G, tau: &G::Tangent) -> G {
        x.lplus(tau)
    }
    fn minus<S: Real, G: LieGroup<S>>(y: &G, x: &G) -> G::Tangent {
        y.lminus(x)
    }
    fn plus_jacobians<S: Real, G: LieGroup<S>>(x: &G, tau: &G::Tangent) -> (G::Jac, G::Jac) {
        x.lplus_jacobians(tau)
    }
    fn minus_jacobians<S: Real, G: LieGroup<S>>(y: &G, x: &G) -> (G::Jac, G::Jac) {
        y.lminus_jacobians(x)
    }
}
