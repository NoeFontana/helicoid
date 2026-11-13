//! Control Jacobians by backward accumulation, `O(K)` per sample (`NUMERICS.md` §21, `0114`).

use helicoid::{Jac, LieGroup, Side, Tangent};
use helicoid_linalg::Real;

use crate::order::{Order, ValidOrder};
use crate::segment::{Kinematics, Segment};
use crate::weights::Weights;

/// The Jacobians of a spline's value, velocity and acceleration with respect to each control.
///
/// Row `m` is with respect to control `m`, perturbed on the evaluation's side. A derivative above
/// the weights' `D` is zero.
#[derive(Clone, Copy)]
pub struct ControlJacobians<S: Real, G: LieGroup<S>, const K: usize> {
    /// `∂X/∂δ_m`, the value perturbed on the same side as the controls.
    pub value: [G::Jac; K],
    /// `∂ω/∂δ_m`.
    pub vel: [G::Jac; K],
    /// `∂ω̇/∂δ_m`.
    pub acc: [G::Jac; K],
}

/// A [`Segment`] plus what every Jacobian sample on it shares: `J_r⁻¹(d_j)`, `J_l⁻¹(d_j)` and
/// `Ad_{X_m⁻¹}` for the left side, formed once per segment.
///
/// Separate from `Segment` so that value-only callers never pay for the inverse Jacobians
/// (`0114` decision 5).
#[derive(Clone, Copy)]
pub struct JacSegment<S: Real, G: LieGroup<S>, const K: usize> {
    seg: Segment<S, G, K>,
    jr_inv: [G::Jac; K],
    jl_inv: [G::Jac; K],
    ad_inv: [G::Jac; K],
}

impl<S: Real, G: LieGroup<S>, const K: usize> JacSegment<S, G, K>
where
    Order<K>: ValidOrder,
{
    /// The segment whose controls are `window`, with its Jacobian factors.
    ///
    /// # Domain
    ///
    /// As [`Segment::new`]: any window.
    #[inline]
    pub fn new(window: &[G; K]) -> Self {
        let seg = Segment::new(window);
        let mut jr_inv = [G::Jac::identity(); K];
        let mut jl_inv = [G::Jac::identity(); K];
        for j in 1..K {
            jr_inv[j] = G::jr_inv(&seg.d[j]);
            jl_inv[j] = G::jl_inv(&seg.d[j]);
        }
        let ad_inv = core::array::from_fn(|m| window[m].inverse().adjoint());
        Self {
            seg,
            jr_inv,
            jl_inv,
            ad_inv,
        }
    }

    /// The underlying segment, for value-only samples on the same window.
    #[inline]
    pub fn segment(&self) -> &Segment<S, G, K> {
        &self.seg
    }

    /// The kinematics and their Jacobians with respect to every control, on side `Sd`.
    ///
    /// `::<Right, _>` perturbs each control as `X_m Exp(δ_m)` and the value as `X Exp(ε)`: one
    /// forward pass keeps `R_j`, `ω^(j)`, `ω̇^(j)`, one backward pass carries `T_j` and `C_j`.
    /// `::<Left, _>` is derived from it, `Ad_X J Ad_{X_m⁻¹}` for the value and `J Ad_{X_m⁻¹}`
    /// for the velocity and acceleration, which are tangents with no output side.
    pub fn kinematics_jacobians<Sd: Side, const D: usize>(
        &self,
        w: &Weights<S, K, D>,
    ) -> (Kinematics<S, G>, ControlJacobians<S, G, K>) {
        let d = &self.seg.d;
        let id = G::Jac::identity();
        let zero_jac = id.scale(S::zero());
        let zero = G::Tangent::zero();

        // Forward: the recursion of `Segment::kinematics`, keeping what the partials read.
        let mut r = [id; K];
        let mut vel_before = [zero; K];
        let mut acc_before = [zero; K];
        let mut vel_at = [zero; K];
        let mut x = self.seg.x0;
        let (mut vel, mut acc) = (zero, zero);
        for j in 1..K {
            let a = G::exp(&d[j].scale(w.lambda[j]));
            x = x * a;
            // `R_j` carries the value's Jacobian too (`T_j`), so it is formed at every `D`.
            r[j] = a.inverse().adjoint();
            if D >= 1 {
                vel_before[j] = vel;
                acc_before[j] = acc;
                vel = r[j].apply(&vel).add(&d[j].scale(w.dot[j]));
                vel_at[j] = vel;
                if D >= 2 {
                    let coriolis = G::ad(&vel).apply(&d[j]).scale(w.dot[j]);
                    acc = r[j].apply(&acc).add(&d[j].scale(w.ddot[j])).add(&coriolis);
                }
            }
        }

        // Backward: `∂/∂d_j` of each output, carried by `T_j` (and `C_j` into `ω̇`).
        let mut dx = [zero_jac; K];
        let mut dv = [zero_jac; K];
        let mut da = [zero_jac; K];
        let (mut t, mut c) = (id, zero_jac);
        for j in (1..K).rev() {
            let jj = G::jr(&d[j].scale(w.lambda[j]));
            dx[j] = t.mul(&jj).scale(w.lambda[j]);
            if D >= 1 {
                let pw = G::ad(&r[j].apply(&vel_before[j]))
                    .mul(&jj)
                    .scale(w.lambda[j])
                    .add(&id.scale(w.dot[j]));
                dv[j] = t.mul(&pw);
                if D >= 2 {
                    let ad_d = G::ad(&d[j]);
                    let pa = G::ad(&r[j].apply(&acc_before[j]))
                        .mul(&jj)
                        .scale(w.lambda[j])
                        .add(&id.scale(w.ddot[j]))
                        .add(&G::ad(&vel_at[j]).scale(w.dot[j]))
                        .add(&ad_d.mul(&pw).scale(-w.dot[j]));
                    da[j] = c.mul(&pw).add(&t.mul(&pa));
                    // `C_{j-1} = C_j R_j − λ̇_j T_j ad(d_j) R_j`, read before `T` moves.
                    c = c.mul(&r[j]).add(&t.mul(&ad_d).mul(&r[j]).scale(-w.dot[j]));
                }
            }
            t = t.mul(&r[j]);
        }

        // Per control: the two increments that read it, `∂d_m/∂δ_m = J_r⁻¹(d_m)` and
        // `∂d_{m+1}/∂δ_m = −J_l⁻¹(d_{m+1})`; the term at `m + 1` added to the term at `m`.
        let mut jac = ControlJacobians {
            value: [zero_jac; K],
            vel: [zero_jac; K],
            acc: [zero_jac; K],
        };
        for m in 0..K {
            let mut terms: [Option<(usize, G::Jac)>; 2] = [None, None];
            if m >= 1 {
                terms[0] = Some((m, self.jr_inv[m]));
            }
            if m + 1 < K {
                terms[1] = Some((m + 1, self.jl_inv[m + 1].neg()));
            }
            let mut value = if m == 0 { t } else { zero_jac };
            let (mut v, mut a) = (zero_jac, zero_jac);
            let mut first = m != 0;
            for (j, dd) in terms.into_iter().flatten() {
                let (tx, tv, ta) = (dx[j].mul(&dd), dv[j].mul(&dd), da[j].mul(&dd));
                if first {
                    // Start from the first term, not `0 + term`, which would round the same but
                    // costs a sum.
                    value = tx;
                    v = tv;
                    a = ta;
                    first = false;
                } else {
                    value = value.add(&tx);
                    v = v.add(&tv);
                    a = a.add(&ta);
                }
            }
            jac.value[m] = value;
            jac.vel[m] = v;
            jac.acc[m] = a;
        }

        if !Sd::IS_RIGHT {
            let ad_x: G::Jac = x.adjoint();
            for m in 0..K {
                jac.value[m] = ad_x.mul(&jac.value[m]).mul(&self.ad_inv[m]);
                let (v, a): (G::Jac, G::Jac) = (jac.vel[m], jac.acc[m]);
                jac.vel[m] = v.mul(&self.ad_inv[m]);
                jac.acc[m] = a.mul(&self.ad_inv[m]);
            }
        }
        (Kinematics { value: x, vel, acc }, jac)
    }
}
