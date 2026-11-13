//! One segment of a uniform cumulative B-spline: value and kinematics (`NUMERICS.md` §21).

use core::marker::PhantomData;

use helicoid::{Jac, LieGroup, Tangent};
use helicoid_linalg::Real;

use crate::order::{Order, ValidOrder};
use crate::weights::Weights;

/// The value of a spline and its right-trivialized derivatives at one `u`.
///
/// `vel` is `(X⁻¹ dX/dt)^∨` and `acc` its time derivative: body rate and body angular acceleration
/// on SO(3), world-frame `ṗ` and `p̈` on Rⁿ, each component on a product. A derivative above the
/// weights' `D` is zero.
#[derive(Clone, Copy, Debug)]
pub struct Kinematics<S: Real, G: LieGroup<S>> {
    /// `X(u)`.
    pub value: G,
    /// `ω`, right-trivialized velocity.
    pub vel: G::Tangent,
    /// `ω̇`, its time derivative.
    pub acc: G::Tangent,
}

/// A segment's `x0` and increments `d_j = Log(X_{j-1}⁻¹ X_j)`, formed once and shared by every
/// sample on it (`0114` decision 5).
#[derive(Clone, Copy, Debug)]
pub struct Segment<S: Real, G: LieGroup<S>, const K: usize> {
    pub(crate) x0: G,
    /// Slot 0 is unused: `K - 1` as an array length needs `generic_const_exprs`.
    pub(crate) d: [G::Tangent; K],
    _s: PhantomData<S>,
}

impl<S: Real, G: LieGroup<S>, const K: usize> Segment<S, G, K>
where
    Order<K>: ValidOrder,
{
    /// The segment whose controls are `window`, in the caller's storage and order.
    ///
    /// # Domain
    ///
    /// Any window: `Log` is total, so near a relative rotation of π the increment follows its
    /// documented behaviour (`NUMERICS.md` §3.2, §12) and there is no failure value.
    #[inline]
    pub fn new(window: &[G; K]) -> Self {
        let mut d = [G::Tangent::zero(); K];
        for j in 1..K {
            d[j] = window[j].rminus(&window[j - 1]);
        }
        Self {
            x0: window[0],
            d,
            _s: PhantomData,
        }
    }

    /// `X(u) = X_0 Exp(λ_1 d_1) ⋯ Exp(λ_{K-1} d_{K-1})`, left to right.
    #[inline]
    pub fn value<const D: usize>(&self, w: &Weights<S, K, D>) -> G {
        let mut x = self.x0;
        for j in 1..K {
            x = x * G::exp(&self.d[j].scale(w.lambda[j]));
        }
        x
    }

    /// The value and its derivatives up to the weights' `D`, by the forward recursion of
    /// `NUMERICS.md` §21. `D = 0` is [`value`](Self::value) with zero derivatives.
    #[inline]
    pub fn kinematics<const D: usize>(&self, w: &Weights<S, K, D>) -> Kinematics<S, G> {
        let zero = G::Tangent::zero();
        let mut x = self.x0;
        let (mut vel, mut acc) = (zero, zero);
        for j in 1..K {
            let a = G::exp(&self.d[j].scale(w.lambda[j]));
            x = x * a;
            if D >= 1 {
                // `R_j = Ad_{A_j⁻¹}`, one adjoint shared by both recursions.
                let r = a.inverse().adjoint();
                vel = r.apply(&vel).add(&self.d[j].scale(w.dot[j]));
                if D >= 2 {
                    let coriolis = G::ad(&vel).apply(&self.d[j]).scale(w.dot[j]);
                    acc = r
                        .apply(&acc)
                        .add(&self.d[j].scale(w.ddot[j]))
                        .add(&coriolis);
                }
            }
        }
        Kinematics { value: x, vel, acc }
    }
}
