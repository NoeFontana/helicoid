//! Which end a perturbation multiplies (`docs/PHASE3.md` §2, `docs/API.md` R1).

use crate::traits::LieGroup;
use helicoid_linalg::Real;

mod sealed {
    pub trait Sealed {}
}

/// The end of a product a perturbation multiplies: [`Right`] (`X · Exp(τ)`, the default) or
/// [`Left`] (`Exp(τ) · X`). Generic code takes `Sd: Side`, never a flag, and each side's Jacobians
/// are expressed in that side's convention (`NUMERICS.md` §2.3).
///
/// **Sealed: [`Right`] and [`Left`] are the only implementations** (`0025`). A group reads its side
/// row by comparing `TypeId`s, which is only sound because this set is closed: an outside impl
/// would get its own `plus`/`minus` and the *`Left`* row of `compose_jacobians`, with no compile
/// error and no runtime signal. Sealing also keeps the choice between `TypeId` and a first-class
/// selector an internal one, which is why `PHASE3.md` §2 can still leave it to the SO(3) PR.
///
/// A downstream implementation does not compile; the delegation below is otherwise complete, so
/// the only error is the unsatisfied supertrait:
///
/// ```compile_fail,E0277
/// use helicoid::{LieGroup, Side};
/// use helicoid_linalg::Real;
///
/// #[derive(Clone, Copy)]
/// struct MyRight;
///
/// impl Side for MyRight {
///     const IS_RIGHT: bool = true;
///     fn plus<S: Real, G: LieGroup<S>>(x: &G, tau: &G::Tangent) -> G {
///         x.rplus(tau)
///     }
///     fn minus<S: Real, G: LieGroup<S>>(y: &G, x: &G) -> G::Tangent {
///         y.rminus(x)
///     }
///     fn plus_jacobians<S: Real, G: LieGroup<S>>(x: &G, tau: &G::Tangent) -> (G::Jac, G::Jac) {
///         x.rplus_jacobians(tau)
///     }
///     fn minus_jacobians<S: Real, G: LieGroup<S>>(y: &G, x: &G) -> (G::Jac, G::Jac) {
///         y.rminus_jacobians(x)
///     }
/// }
/// ```
pub trait Side: sealed::Sealed + Copy + 'static {
    /// `true` for [`Right`], `false` for [`Left`]: how a group reads its side row of
    /// `NUMERICS.md` §2.3 where the two differ (`SO3::compose_jacobians`).
    ///
    /// This is the "first-class selector" the trait docs weigh against a `TypeId` comparison, and
    /// the SO(3) PR took it (`PHASE3.md` §2 left the choice here). It is an associated **const**,
    /// so `match Sd::IS_RIGHT` is resolved at monomorphization and neither branch survives into
    /// the emitted code — a `TypeId` compare is a runtime call on a 16-byte value that LLVM folds
    /// only after inlining, and it needs the `'static` bound to stay. A group that reads the same
    /// Jacobian on both sides, such as `Rn`, never mentions it.
    const IS_RIGHT: bool;
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

impl sealed::Sealed for Right {}
impl sealed::Sealed for Left {}

impl Side for Right {
    const IS_RIGHT: bool = true;
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
    const IS_RIGHT: bool = false;
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
