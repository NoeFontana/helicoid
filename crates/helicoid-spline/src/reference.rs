//! Reference twins (D6, `NUMERICS.md` §14): the same quantities by an independent route, `O(K²)`,
//! for the `*_matches_reference` proptests. Never on a hot path.

use helicoid::{Jac, LieGroup, Tangent};
use helicoid_linalg::Real;

use crate::order::{Order, ValidOrder};
use crate::segment::Kinematics;
use crate::weights::Weights;

fn factors<S: Real, G: LieGroup<S>, const K: usize, const D: usize>(
    window: &[G; K],
    w: &Weights<S, K, D>,
) -> ([G::Tangent; K], [G; K]) {
    let mut d = [G::Tangent::zero(); K];
    let mut a = [G::identity(); K];
    for j in 1..K {
        // `Log(X_{j-1}⁻¹ X_j)` spelled out, not `rminus`.
        d[j] = (window[j - 1].inverse() * window[j]).log();
        a[j] = G::exp(&d[j].scale(w.lambda[j]));
    }
    (d, a)
}

/// `P = A_{from} ⋯ A_{to}`, empty products the identity.
fn product<S: Real, G: LieGroup<S>, const K: usize>(a: &[G; K], from: usize, to: usize) -> G {
    let mut p = G::identity();
    for f in a.iter().take(to + 1).skip(from) {
        p = p * *f;
    }
    p
}

/// `X_0 (A_1 (A_2 ⋯ A_{K-1}))`: the value, associated right to left.
pub fn spline_value<S: Real, G: LieGroup<S>, const K: usize, const D: usize>(
    window: &[G; K],
    w: &Weights<S, K, D>,
) -> G
where
    Order<K>: ValidOrder,
{
    let (_, a) = factors(window, w);
    let mut p = G::identity();
    for j in (1..K).rev() {
        p = a[j] * p;
    }
    window[0] * p
}

/// The kinematics as closed-form sums over explicit suffix products, not the recursion:
/// `ω = Σ_j Ad_{P_j⁻¹} λ̇_j d_j` and `ω̇ = Σ_j Ad_{P_j⁻¹} (λ̈_j d_j + λ̇_j ad(ω^(j)) d_j)`, with
/// `ω^(j)` the partial sum up to `j` and `P_j = A_{j+1} ⋯ A_{K-1}` (`docs/maths/splines.md` BS.6).
pub fn spline_kinematics<S: Real, G: LieGroup<S>, const K: usize, const D: usize>(
    window: &[G; K],
    w: &Weights<S, K, D>,
) -> Kinematics<S, G>
where
    Order<K>: ValidOrder,
{
    let (d, a) = factors(window, w);
    // `ω^(j) = Σ_{i ≤ j} Ad_{(A_{i+1} ⋯ A_j)⁻¹} λ̇_i d_i`.
    let partial = |j: usize| {
        let mut s = G::Tangent::zero();
        for (i, (di, doti)) in d.iter().zip(w.dot).enumerate().take(j + 1).skip(1) {
            let ad = product(&a, i + 1, j).inverse().adjoint();
            s = s.add(&ad.apply(&di.scale(doti)));
        }
        s
    };
    let (mut vel, mut acc) = (G::Tangent::zero(), G::Tangent::zero());
    for (j, dj) in d.iter().enumerate().skip(1) {
        let t = product(&a, j + 1, K - 1).inverse().adjoint();
        vel = vel.add(&t.apply(&dj.scale(w.dot[j])));
        let own = dj
            .scale(w.ddot[j])
            .add(&G::ad(&partial(j)).apply(dj).scale(w.dot[j]));
        acc = acc.add(&t.apply(&own));
    }
    Kinematics {
        value: spline_value(window, w),
        vel,
        acc,
    }
}

/// The control Jacobians by **forward** accumulation, one pass per control, `O(K²)`: the chain
/// rule pushed through the recursion from `X_0` up, where the kernel pulls it back from the output
/// (`docs/maths/splines.md` BS.6). Same side conversion as the kernel.
pub fn spline_jacobians<
    S: Real,
    G: LieGroup<S>,
    Sd: helicoid::Side,
    const K: usize,
    const D: usize,
>(
    window: &[G; K],
    w: &Weights<S, K, D>,
) -> crate::ControlJacobians<S, G, K>
where
    Order<K>: ValidOrder,
{
    let (d, a) = factors(window, w);
    let id = G::Jac::identity();
    let zero = id.scale(S::zero());
    let mut out = crate::ControlJacobians {
        value: [zero; K],
        vel: [zero; K],
        acc: [zero; K],
    };
    for m in 0..K {
        // `∂d_j/∂δ_m`: `J_r⁻¹(d_m)` at `j = m`, `−J_l⁻¹(d_{m+1})` at `j = m + 1`, else 0.
        let dd = |j: usize| {
            if j == m {
                G::jr_inv(&d[j])
            } else if j == m + 1 {
                G::jl_inv(&d[j]).neg()
            } else {
                zero
            }
        };
        let mut gx = if m == 0 { id } else { zero };
        let (mut gv, mut ga) = (zero, zero);
        let (mut vel, mut acc) = (G::Tangent::zero(), G::Tangent::zero());
        for j in 1..K {
            let r = a[j].inverse().adjoint();
            let jj = G::jr(&d[j].scale(w.lambda[j]));
            let ddj = dd(j);
            let vel_next = r.apply(&vel).add(&d[j].scale(w.dot[j]));
            let pw = G::ad(&r.apply(&vel))
                .mul(&jj)
                .scale(w.lambda[j])
                .add(&id.scale(w.dot[j]));
            let pa = G::ad(&r.apply(&acc))
                .mul(&jj)
                .scale(w.lambda[j])
                .add(&id.scale(w.ddot[j]))
                .add(&G::ad(&vel_next).scale(w.dot[j]))
                .add(&G::ad(&d[j]).mul(&pw).scale(-w.dot[j]));
            let ga_next = r
                .mul(&ga)
                .add(&G::ad(&d[j]).mul(&r).mul(&gv).scale(-w.dot[j]))
                .add(&pa.mul(&ddj));
            gv = r.mul(&gv).add(&pw.mul(&ddj));
            ga = ga_next;
            gx = r.mul(&gx).add(&jj.mul(&ddj).scale(w.lambda[j]));
            let coriolis = G::ad(&vel_next).apply(&d[j]).scale(w.dot[j]);
            acc = r.apply(&acc).add(&d[j].scale(w.ddot[j])).add(&coriolis);
            vel = vel_next;
        }
        out.value[m] = gx;
        out.vel[m] = if D >= 1 { gv } else { zero };
        out.acc[m] = if D >= 2 { ga } else { zero };
    }
    if !Sd::IS_RIGHT {
        let ad_x: G::Jac = spline_value(window, w).adjoint();
        for (m, xm) in window.iter().enumerate() {
            let ad_inv = xm.inverse().adjoint();
            out.value[m] = ad_x.mul(&out.value[m]).mul(&ad_inv);
            let (v, a): (G::Jac, G::Jac) = (out.vel[m], out.acc[m]);
            out.vel[m] = v.mul(&ad_inv);
            out.acc[m] = a.mul(&ad_inv);
        }
    }
    out
}
