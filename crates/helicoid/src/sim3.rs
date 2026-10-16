//! Sim(3): rotation, translation and scale (`docs/PHASE5.md` §3, `NUMERICS.md` §9, `0070`).
//!
//! Stored `(q, t, σ)` with scale `s = e^σ`, so the inverse negates `σ` exactly. Tangents are
//! `[φ; ρ; σ]`, rotation first. The Jacobians are [`Sim3Jac`], the structured
//! `[A 0 0; B C f; 0 0 e]` every Sim(3) Jacobian has (`docs/maths/sim3.md` SM.9), closed under
//! product and inverse; the coefficients are `coeffs::sim3`'s group.

use crate::coeffs::{
    b_takes_exact_arm, sim3_coeffs, sim3_coeffs_scaled, sim3_coeffs_scaled_pair, sim3_v,
    sim3_v_scaled, unscale, Sim3Coeffs,
};
use crate::quat::Quat;
use crate::side::Side;
use crate::so3::{hat_mul, mul_hat, norm_sq, SO3Tangent, SO3};
use crate::traits::{tie_dof, Jac, LieGroup, Tangent};
use core::array;
use core::ops::Mul;
use helicoid_linalg::{
    hat, Blend, Mask, Mat3, Matrix, Point, Point3, Real, StridedMut, Vec3, Vector,
};

/// An element of Sim(3): `X p = s R p + t`, `s = e^σ` (`NUMERICS.md` §9).
#[derive(Clone, Copy, Debug)]
pub struct Sim3<S> {
    q: Quat<S>,
    t: Vec3<S>,
    sigma: S,
}

/// A tangent vector of [`Sim3`]: `[φ; ρ; σ]`, rotation first (`0002`).
#[derive(Clone, Copy, Debug)]
pub struct Sim3Tangent<S> {
    /// The rotation block.
    pub phi: Vec3<S>,
    /// The translation block.
    pub rho: Vec3<S>,
    /// The log-scale.
    pub sigma: S,
}

/// A Sim(3) Jacobian, `[A 0 0; B C f; 0 0 e]` on `[φ; ρ; σ]` (`0070` decision 1): `A`, `B`, `C`
/// `3 x 3`, `f` a column, `e` a scalar. `Ad`, `ad`, `J_l`, `J_r`, their inverses and every row of
/// `NUMERICS.md` §2.3 built from them have it, and it is closed under [`mul`](Jac::mul) (`117`
/// multiplications against `343` dense) and [`inverse`](Jac::inverse).
#[derive(Clone, Copy, Debug)]
pub struct Sim3Jac<S> {
    pub(crate) a: Mat3<S>,
    pub(crate) b: Mat3<S>,
    pub(crate) c: Mat3<S>,
    pub(crate) f: Vec3<S>,
    pub(crate) e: S,
}

impl<S: Real> Blend<S> for Sim3<S> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self {
            q: Quat::blend(m, t.q, f.q),
            t: Vector::blend(m, t.t, f.t),
            sigma: S::select(m, t.sigma, f.sigma),
        }
    }
}

impl<S: Real> Blend<S> for Sim3Tangent<S> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self {
            phi: Vector::blend(m, t.phi, f.phi),
            rho: Vector::blend(m, t.rho, f.rho),
            sigma: S::select(m, t.sigma, f.sigma),
        }
    }
}

impl<S: Real> Blend<S> for Sim3Jac<S> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self {
            a: Mat3::blend(m, t.a, f.a),
            b: Mat3::blend(m, t.b, f.b),
            c: Mat3::blend(m, t.c, f.c),
            f: Vector::blend(m, t.f, f.f),
            e: S::select(m, t.e, f.e),
        }
    }
}

impl<S: Real> Sim3Tangent<S> {
    #[inline]
    fn comps(&self) -> [S; 7] {
        let ([a, b, c], [d, e, f]) = (self.phi.0, self.rho.0);
        [a, b, c, d, e, f, self.sigma]
    }

    #[inline]
    fn from_comps(v: [S; 7]) -> Self {
        let [a, b, c, d, e, f, g] = v;
        Self {
            phi: Vector([a, b, c]),
            rho: Vector([d, e, f]),
            sigma: g,
        }
    }
}

impl<S: Real> Tangent<S> for Sim3Tangent<S> {
    const DOF: usize = 7;
    #[inline]
    fn zero() -> Self {
        Self::from_comps([S::zero(); 7])
    }
    #[inline]
    fn add(&self, o: &Self) -> Self {
        let (a, b) = (self.comps(), o.comps());
        Self::from_comps(array::from_fn(|i| a[i] + b[i]))
    }
    #[inline]
    fn sub(&self, o: &Self) -> Self {
        let (a, b) = (self.comps(), o.comps());
        Self::from_comps(array::from_fn(|i| a[i] - b[i]))
    }
    #[inline]
    fn neg(&self) -> Self {
        Self::from_comps(self.comps().map(|x| -x))
    }
    #[inline]
    fn scale(&self, k: S) -> Self {
        Self::from_comps(self.comps().map(|x| k * x))
    }
    #[inline]
    fn dot_acc(&self, o: &Self, acc: S) -> S {
        let (a, b) = (self.comps(), o.comps());
        (0..7).fold(acc, |s, i| s + a[i] * b[i])
    }
    #[inline]
    fn write_dense(&self, out: &mut [S]) {
        debug_assert!(out.len() == 7, "Tangent::write_dense: wrong length");
        for (o, v) in out.iter_mut().zip(self.comps()) {
            *o = v;
        }
    }
    #[inline]
    fn read_dense(src: &[S]) -> Self {
        debug_assert!(src.len() == 7, "Tangent::read_dense: wrong length");
        // A short `src` reads NaN, not `+0`, for `RnTangent::read_dense`'s reason (D11).
        Self::from_comps(array::from_fn(|i| {
            src.get(i).copied().unwrap_or_else(|| S::zero() / S::zero())
        }))
    }
}

impl<S: Real> Sim3<S> {
    /// `(R, t, σ)`, the rotation moved in.
    #[inline]
    pub fn from_parts(r: SO3<S>, t: Vec3<S>, sigma: S) -> Self {
        Self {
            q: r.quat(),
            t,
            sigma,
        }
    }

    /// `(R, t, σ)`.
    #[inline]
    pub fn parts(&self) -> (SO3<S>, Vec3<S>, S) {
        (self.rotation(), self.t, self.sigma)
    }

    /// The rotation.
    #[inline]
    pub fn rotation(&self) -> SO3<S> {
        SO3::from_quat_unchecked(self.q)
    }

    /// The translation.
    #[inline]
    pub fn translation(&self) -> Vec3<S> {
        self.t
    }

    /// `σ`, the stored log-scale.
    #[inline]
    pub fn log_scale(&self) -> S {
        self.sigma
    }

    /// `s = e^σ`.
    ///
    /// # Domain
    ///
    /// `σ < ln MAX` for a finite scale (`NUMERICS.md` §12).
    #[inline]
    pub fn scale(&self) -> S {
        self.sigma.exp()
    }

    /// The quaternion's Newton step, in place (`API.md` R2): only the rotation drifts, since `t`
    /// and `σ` carry no constraint.
    #[inline]
    pub fn renormalize(&mut self) {
        self.q.renormalize();
    }

    /// `self · other⁻¹ = (R R_o⁻¹, t - e^{σ - σ_o} R R_o⁻¹ t_o, σ - σ_o)`, the relative rotation
    /// formed once and the scales subtracted before the exponential.
    #[inline]
    pub fn mul_inv(&self, other: &Self) -> Self {
        let q = self.q * other.q.conjugate();
        let ds = self.sigma - other.sigma;
        let rt = SO3::from_quat_unchecked(q).act(other.t);
        Self {
            q,
            t: self.t - rt.scale(ds.exp()),
            sigma: ds,
        }
    }

    /// `other⁻¹ · self = (R_o⁻¹ R, e^{-σ_o} R_o⁻¹ (t - t_o), σ - σ_o)`, subtracting before the
    /// rotation, where it is exact for nearby elements.
    #[inline]
    pub fn inv_mul(&self, other: &Self) -> Self {
        let oc = other.q.conjugate();
        let d = SO3::from_quat_unchecked(oc).act(self.t - other.t);
        Self {
            q: oc * self.q,
            t: d.scale((-other.sigma).exp()),
            sigma: self.sigma - other.sigma,
        }
    }

    /// `X p = s R p + t`; public as `Mul<Point3>`, the one spelling (`API.md` R1, R3).
    #[inline]
    pub(crate) fn act(&self, p: Point3<S>) -> Point3<S> {
        let rp = self.rotation().act(Vector(p.0));
        Point((rp.scale(self.scale()) + self.t).0)
    }

    /// `X p` over a slice, the rotation matrix and the scale formed once.
    #[inline]
    pub fn act_many(&self, pts: &mut [Point3<S>]) {
        let sr = self.rotation().to_matrix().scale(self.scale());
        for p in pts.iter_mut() {
            *p = Point((sr * Vector(p.0) + self.t).0);
        }
    }

    /// `(∂(X p)/∂X, ∂(X p)/∂p)` in the convention of `Sd` (`docs/maths/sim3.md` SM.12): right
    /// `[-sR[p]_×, sR, sRp]`, left `[-[X p]_×, I, X p]`, and `sR`. A dense `3 x 7` for `SEn3`'s
    /// reason (`0005`).
    #[inline]
    pub fn act_jacobians<Sd: Side>(&self, p: Point3<S>) -> (Matrix<S, 3, 7>, Mat3<S>) {
        let sr = self.rotation().to_matrix().scale(self.scale());
        let (rot, trans, scale) = match Sd::IS_RIGHT {
            true => (-mul_hat(&sr, Vector(p.0)), sr, sr * Vector(p.0)),
            false => {
                let xp = Vector(self.act(p).0);
                (-hat(xp), Mat3::identity(), xp)
            }
        };
        let cols = array::from_fn(|c| match c {
            0..=2 => rot.col(c),
            3..=5 => trans.col(c - 3),
            _ => scale,
        });
        (Matrix::from_cols(cols), sr)
    }
}

impl<S: Real> Mul for Sim3<S> {
    type Output = Self;
    /// `XY = (R_X R_Y, s_X R_X t_Y + t_X, σ_X + σ_Y)` (SM.2).
    #[inline]
    fn mul(self, o: Self) -> Self {
        let rt = self.rotation().act(o.t);
        Self {
            q: self.q * o.q,
            t: rt.scale(self.scale()) + self.t,
            sigma: self.sigma + o.sigma,
        }
    }
}

impl<S: Real> Mul<Point3<S>> for Sim3<S> {
    type Output = Point3<S>;
    #[inline]
    fn mul(self, p: Point3<S>) -> Point3<S> {
        self.act(p)
    }
}

/// `u vᵀ`.
#[inline]
fn outer<S: Real>(u: Vec3<S>, v: Vec3<S>) -> Mat3<S> {
    Matrix::from_cols(v.0.map(|x| u.scale(x)))
}

/// `x₀ I + x₁ W + x₂ W²` applied to `v` by two cross products, `W = [φ]_×`.
#[inline]
fn poly_apply<S: Real>(phi: Vec3<S>, x: [S; 3], v: Vec3<S>) -> Vec3<S> {
    let mu = phi.cross(v);
    (v.scale(x[0]) + mu.scale(x[1])) + phi.cross(mu).scale(x[2])
}

/// `Vρ`, `V = v₀I + v₁W + v₂W²` at `z = ‖φ‖²`: by two cross products, or where `b` is on its exact
/// arm `χρ + v₁ (φ × ρ) + v₂ (φ·ρ) φ` with `χ = v₀ − zv₂` from the kernel's `χΔ`, formed without that
/// difference, which cancels to `sin θ/θ` at `σ = 0` near `π`. `SEn3::exp`'s far arm
/// (`0095`, draft; `0096`, draft), one branch on the same mask; `Δ`'s division is in the arm that
/// reads it, which measured 1 ns cheaper on the far bands than in the kernel.
#[inline]
fn v_apply<S: Real>(
    phi: Vec3<S>,
    z: S,
    sigma: S,
    v: [S; 3],
    chi_delta: S,
    rho: Vec3<S>,
) -> Vec3<S> {
    S::branch(
        b_takes_exact_arm(z),
        move || {
            let chi = chi_delta / (sigma * sigma + z);
            (rho.scale(chi) + phi.cross(rho).scale(v[1])) + phi.scale(v[2] * phi.dot(rho))
        },
        move || poly_apply(phi, v, rho),
    )
}

/// `x₀ I + x₁ W + x₂ W²` as a matrix.
#[inline]
fn poly<S: Real>(phi: Vec3<S>, x: [S; 3]) -> Mat3<S> {
    let w = hat(phi);
    Mat3::identity().scale(x[0]) + w.scale(x[1]) + mul_hat(&w, phi).scale(x[2])
}

/// `[î₀, î₁, î₂]` of `V̂⁻¹` from `V̂`'s `v̂` (SM.8): `χ = v₀ - z v₂`, `Λ = χ² + z v₁²`.
#[inline]
fn inv_coeffs<S: Real>(v: [S; 3], z: S) -> [S; 3] {
    let [v0, v1, v2] = v;
    let chi = v0 - z * v2;
    let lambda = chi * chi + z * (v1 * v1);
    [
        S::one() / v0,
        -(v1 / lambda),
        (v1 * v1 - v2 * chi) / (v0 * lambda),
    ]
}

/// `J_l(τ)` from the group's coefficients `c` at `(σ, ‖φ‖²)`, unscaled (SM.10):
/// `[J_l(φ) 0 0; Q V −Γ₂(M)ρ; 0 0 1]`, `y = Vρ`,
/// `Q = [y]_× J_l(φ) + 2(∂_z v₁ μ + ∂_z v₂ ν)φᵀ − v₁Ξ + v₂(ΞW − 2WΞ)`, `μ = φ × ρ`, `ν = φ × μ`.
///
/// SO(3)'s `J_l(φ)` is formed here, not passed in: as an argument `jr` paid 2–9% to carry it
/// through memory into the out-of-line call (`0073`).
#[inline]
fn jl_blocks<S: Real>(phi: Vec3<S>, rho: Vec3<S>, c: &Sim3Coeffs<S>) -> Sim3Jac<S> {
    jl_blocks_with(phi, rho, c, SO3::jl(&SO3Tangent { phi }))
}

/// [`jl_blocks`] given SO(3)'s `J_l(φ)`.
#[inline(always)]
fn jl_blocks_with<S: Real>(
    phi: Vec3<S>,
    rho: Vec3<S>,
    c: &Sim3Coeffs<S>,
    jl_phi: Mat3<S>,
) -> Sim3Jac<S> {
    let (w, xi) = (hat(phi), hat(rho));
    let mu = phi.cross(rho);
    let nu = phi.cross(mu);
    let y = (rho.scale(c.v[0]) + mu.scale(c.v[1])) + nu.scale(c.v[2]);
    let gamma = (rho.scale(c.w[0]) + mu.scale(c.w[1])) + nu.scale(c.w[2]);
    let two = S::lit(2.0);
    let dy = mu.scale(two * c.dv[0]) + nu.scale(two * c.dv[1]);
    let q = ((hat_mul(y, &jl_phi) + outer(dy, phi)) - xi.scale(c.v[1]))
        + (xi * w - (w * xi).scale(two)).scale(c.v[2]);
    Sim3Jac {
        a: jl_phi,
        b: q,
        c: poly(phi, c.v),
        f: -gamma,
        e: S::one(),
    }
}

/// SO(3)'s `(J_l(φ), J_l⁻¹(φ))`. At `-φ` they are these transposed, to the bit
/// (`so3_tests::jl_is_jr_transposed_to_the_bit`).
#[inline]
fn so3_jl_pair<S: Real>(phi: Vec3<S>) -> (Mat3<S>, Mat3<S>) {
    let t = SO3Tangent { phi };
    (SO3::jl(&t), SO3::jl_inv(&t))
}

/// `(J_l(τ), J_l(τ)⁻¹)` by blocks (SM.10(c)): `A⁻¹` from SO(3)'s closed form, `C⁻¹ = V⁻¹` from
/// SM.8 on the scaled `v̂` and `e^{-σ}` applied last, so `Λ` does not overflow before `V⁻¹` does.
#[inline]
fn jl_and_inverse<S: Real>(phi: Vec3<S>, rho: Vec3<S>, sigma: S) -> (Sim3Jac<S>, Sim3Jac<S>) {
    let z = norm_sq(phi);
    jl_and_inverse_from(
        phi,
        rho,
        sigma,
        z,
        sim3_coeffs_scaled(sigma, z),
        so3_jl_pair(phi),
    )
}

/// [`jl_and_inverse`] from the group's `(ĉ, scaled)` at `(σ, z = ‖φ‖²)` and [`so3_jl_pair`]`(φ)`.
#[inline]
fn jl_and_inverse_from<S: Real>(
    phi: Vec3<S>,
    rho: Vec3<S>,
    sigma: S,
    z: S,
    (hat_c, scaled): (Sim3Coeffs<S>, S::Mask),
    (jl_phi, a_inv): (Mat3<S>, Mat3<S>),
) -> (Sim3Jac<S>, Sim3Jac<S>) {
    let j = jl_blocks_with(phi, rho, &unscale(hat_c, sigma, scaled), jl_phi);
    let inv_scale = S::branch(scaled, || (-sigma).exp(), S::one);
    let c_inv = poly(phi, inv_coeffs(hat_c.v, z)).scale(inv_scale);
    let inv = Sim3Jac {
        a: a_inv,
        b: -(c_inv * j.b * a_inv),
        c: c_inv,
        f: -(c_inv * j.f),
        e: S::one(),
    };
    (j, inv)
}

impl<S: Real> Sim3<S> {
    /// `Log`'s `ρ = V⁻¹t` from `v̂` at `(σ, z = ‖φ‖²)` (SM.8), applied as `e^{-σ/2}(V̂⁻¹(e^{-σ/2} t))`
    /// where scaled (SM.7(g)): neither `Λ` nor `V̂⁻¹ t` overflows before the answer does.
    #[inline]
    fn log_rho(&self, phi: Vec3<S>, z: S, (v, scaled): ([S; 3], S::Mask)) -> Vec3<S> {
        let h = S::branch(scaled, || (S::lit(-0.5) * self.sigma).exp(), S::one);
        poly_apply(phi, inv_coeffs(v, z), self.t.scale(h)).scale(h)
    }
}

/// `(Exp τ, J_l(τ))` from one evaluation of the group: `Exp`'s `v` and `χΔ` are the full group's,
/// the same bits as `sim3_v`'s.
#[inline]
fn exp_and_jl<S: Real>(tau: &Sim3Tangent<S>) -> (Sim3<S>, Sim3Jac<S>) {
    let (phi, rho) = (tau.phi, tau.rho);
    let t = SO3Tangent { phi };
    let z = norm_sq(phi);
    let c = sim3_coeffs(tau.sigma, z);
    let x = Sim3 {
        q: SO3::exp(&t).quat(),
        t: v_apply(phi, z, tau.sigma, c.v, c.chi_delta, rho),
        sigma: tau.sigma,
    };
    (x, jl_blocks(phi, rho, &c))
}

/// `(J_l⁻¹(τ), J_r⁻¹(τ))` at `τ = Log x`, each `LieGroup::jl_inv`'s and `jr_inv`'s to the bit
/// (`0073`): the group's coefficients at `σ` and `-σ` from one [`sim3_coeffs_scaled_pair`], whose
/// `v` at `σ` is `Log`'s, and SO(3)'s pair at `φ` transposed for `-φ`. Against `Log` and the two
/// calls this drops two of three coefficient evaluations and one of two SO(3) pairs.
#[inline]
fn log_inverses<S: Real>(x: &Sim3<S>) -> (Sim3Jac<S>, Sim3Jac<S>) {
    let phi = x.rotation().log().phi;
    let z = norm_sq(phi);
    let [pos, neg] = sim3_coeffs_scaled_pair(x.sigma, z);
    let rho = x.log_rho(phi, z, (pos.0.v, pos.1));
    let (jl, jl_inv) = so3_jl_pair(phi);
    let left = jl_and_inverse_from(phi, rho, x.sigma, z, pos, (jl, jl_inv)).1;
    let right = jl_and_inverse_from(
        -phi,
        -rho,
        -x.sigma,
        z,
        neg,
        (jl.transpose(), jl_inv.transpose()),
    )
    .1;
    (left, right)
}

impl<S: Real> LieGroup<S> for Sim3<S> {
    type Tangent = Sim3Tangent<S>;
    type Jac = Sim3Jac<S>;
    const DOF: usize = 7;

    #[inline]
    fn identity() -> Self {
        const { tie_dof::<S, Self>() };
        Self::from_parts(SO3::identity(), Vector([S::zero(); 3]), S::zero())
    }

    /// `X⁻¹ = (R⁻¹, −e^{-σ} R⁻¹ t, −σ)` (SM.2).
    #[inline]
    fn inverse(&self) -> Self {
        let r = self.rotation().inverse();
        Self {
            q: r.quat(),
            t: -r.act(self.t).scale((-self.sigma).exp()),
            sigma: -self.sigma,
        }
    }

    /// `Exp τ = (Exp φ, Vρ, σ)`, `V = v₀I + v₁W + v₂W²` (SM.4): by two cross products, or where
    /// `b` is on its exact arm by `χ = v₀ − zv₂` formed without that cancellation (`0096`, draft).
    #[inline]
    fn exp(tau: &Sim3Tangent<S>) -> Self {
        let phi = tau.phi;
        let z = norm_sq(phi);
        let (v, chi_delta) = sim3_v(tau.sigma, z);
        let q = SO3::exp(&SO3Tangent { phi }).quat();
        Self {
            q,
            t: v_apply(phi, z, tau.sigma, v, chi_delta, tau.rho),
            sigma: tau.sigma,
        }
    }

    /// `Log X = (φ, V⁻¹t, σ)`, `φ` through the quaternion's `atan2` (§3.2), `V⁻¹` by SM.8 on the
    /// scaled `v̂` and applied as `e^{-σ/2}(V̂⁻¹(e^{-σ/2} t))` where scaled (SM.7(g)): neither
    /// `Λ` nor `V̂⁻¹ t` overflows before the answer does.
    ///
    /// # Domain
    ///
    /// `θ ≤ π`, `|σ| < ln MAX` (`NUMERICS.md` §12).
    #[inline]
    fn log(&self) -> Sim3Tangent<S> {
        let phi = self.rotation().log().phi;
        let z = norm_sq(phi);
        Sim3Tangent {
            phi,
            rho: self.log_rho(phi, z, sim3_v_scaled(self.sigma, z)),
            sigma: self.sigma,
        }
    }

    /// `Log(other⁻¹ · self)` through [`inv_mul`](Sim3::inv_mul), `SEn3::rminus`'s reasons.
    #[inline]
    fn rminus(&self, base: &Self) -> Sim3Tangent<S> {
        self.inv_mul(base).log()
    }

    /// `Log(self · other⁻¹)` through [`mul_inv`](Sim3::mul_inv).
    #[inline]
    fn lminus(&self, base: &Self) -> Sim3Tangent<S> {
        self.mul_inv(base).log()
    }

    /// `Ad_X = [R 0 0; [t]_× R  sR  −t; 0 0 1]` (SM.9).
    #[inline]
    fn adjoint(&self) -> Sim3Jac<S> {
        let r = self.rotation().to_matrix();
        Sim3Jac {
            a: r,
            b: hat_mul(self.t, &r),
            c: r.scale(self.scale()),
            f: -self.t,
            e: S::one(),
        }
    }

    /// `ad_τ = [W 0 0; Ξ  W + σI  −ρ; 0 0 0]` (SM.9).
    #[inline]
    fn ad(tau: &Sim3Tangent<S>) -> Sim3Jac<S> {
        let w = hat(tau.phi);
        Sim3Jac {
            a: w,
            b: hat(tau.rho),
            c: w + Mat3::identity().scale(tau.sigma),
            f: -tau.rho,
            e: S::zero(),
        }
    }

    /// `J_r(τ) = J_l(−τ)` (SM.10); `jl` keeps the provided `J_r(−τ)`, which is `J_l(τ)` by the same
    /// one path, as `SEn3::jr`'s note says of its own.
    #[inline]
    fn jr(tau: &Sim3Tangent<S>) -> Sim3Jac<S> {
        let n = tau.neg();
        jl_blocks(n.phi, n.rho, &sim3_coeffs(n.sigma, norm_sq(n.phi)))
    }

    /// `J_r⁻¹(τ) = J_l(−τ)⁻¹` by blocks (SM.10(c)), not by [`Jac::inverse`]: SO(3)'s closed form
    /// for `A⁻¹` and SM.8's for `V⁻¹`, with no `det` to divide by.
    ///
    /// # Domain
    ///
    /// `θ < 2π`, `|σ| < ln MAX` (`NUMERICS.md` §12).
    #[inline]
    fn jr_inv(tau: &Sim3Tangent<S>) -> Sim3Jac<S> {
        let n = tau.neg();
        jl_and_inverse(n.phi, n.rho, n.sigma).1
    }

    /// `(Ad_Exp(τ)⁻¹, J_r(τ))`, `Ad⁻¹` as `Ad_Exp(−τ)` (`NUMERICS.md` §2.3), from one
    /// coefficient evaluation: `Exp`'s `V` is the group's `v`, the same bits (`0073`).
    #[inline]
    fn rplus_jacobians(&self, tau: &Sim3Tangent<S>) -> (Sim3Jac<S>, Sim3Jac<S>) {
        let (x, jr) = exp_and_jl(&tau.neg());
        (x.adjoint(), jr)
    }

    /// `(Ad_Exp(τ), J_l(τ))` (`NUMERICS.md` §2.3), from one coefficient evaluation.
    #[inline]
    fn lplus_jacobians(&self, tau: &Sim3Tangent<S>) -> (Sim3Jac<S>, Sim3Jac<S>) {
        let (x, jl) = exp_and_jl(tau);
        (x.adjoint(), jl)
    }

    /// `(J_r⁻¹(τ), −J_l⁻¹(τ))` at `τ = self ⊖_R base` (`NUMERICS.md` §2.3): `Log`'s coefficients
    /// and both inverses' from one evaluation at `±σ`, each `jr_inv`'s and `jl_inv`'s to the bit
    /// (`0073`).
    #[inline]
    fn rminus_jacobians(&self, base: &Self) -> (Sim3Jac<S>, Sim3Jac<S>) {
        let (jl_inv, jr_inv) = log_inverses(&self.inv_mul(base));
        (jr_inv, jl_inv.neg())
    }

    /// `(J_l⁻¹(τ), −J_r⁻¹(τ))` at `τ = self ⊖_L base` (`NUMERICS.md` §2.3), as
    /// [`rminus_jacobians`](LieGroup::rminus_jacobians).
    #[inline]
    fn lminus_jacobians(&self, base: &Self) -> (Sim3Jac<S>, Sim3Jac<S>) {
        let (jl_inv, jr_inv) = log_inverses(&self.mul_inv(base));
        (jl_inv, jr_inv.neg())
    }

    /// Right `(Ad_Y⁻¹, I)`, left `(I, Ad_X)` (`NUMERICS.md` §2.3).
    #[inline]
    fn compose_jacobians<Sd: Side>(&self, rhs: &Self) -> (Sim3Jac<S>, Sim3Jac<S>) {
        match Sd::IS_RIGHT {
            true => (rhs.inverse().adjoint(), Sim3Jac::identity()),
            false => (Sim3Jac::identity(), self.adjoint()),
        }
    }

    /// Right `−Ad_X`, left `−Ad_X⁻¹` (`NUMERICS.md` §2.3).
    #[inline]
    fn inverse_jacobian<Sd: Side>(&self) -> Sim3Jac<S> {
        match Sd::IS_RIGHT {
            true => self.adjoint().neg(),
            false => self.inverse().adjoint().neg(),
        }
    }
}

/// `J` applied to a dense `[φ; ρ; σ]` held in `v[0..7]`, the structure's own product.
#[inline]
fn apply7<S: Real, const D: usize>(j: &Sim3Jac<S>, v: &Vector<S, D>) -> Vector<S, D> {
    let p = Vector([v.0[0], v.0[1], v.0[2]]);
    let r = Vector([v.0[3], v.0[4], v.0[5]]);
    let s = v.0[6];
    let top = j.a * p;
    let mid = (j.b * p + j.c * r) + j.f.scale(s);
    Vector(array::from_fn(|i| match i {
        0..=2 => top.0[i],
        3..=5 => mid.0[i - 3],
        _ => j.e * s,
    }))
}

impl<S: Real> Jac<S, Sim3Tangent<S>> for Sim3Jac<S> {
    #[inline]
    fn identity() -> Self {
        let zero = Matrix::from_cols([Vector([S::zero(); 3]); 3]);
        Self {
            a: Mat3::identity(),
            b: zero,
            c: Mat3::identity(),
            f: Vector([S::zero(); 3]),
            e: S::one(),
        }
    }
    /// `(AA', BA' + CB', CC', Cf' + f e', e e')`.
    #[inline]
    fn mul(&self, o: &Self) -> Self {
        Self {
            a: self.a * o.a,
            b: self.b * o.a + self.c * o.b,
            c: self.c * o.c,
            f: self.c * o.f + self.f.scale(o.e),
            e: self.e * o.e,
        }
    }
    /// `(A⁻¹, −C⁻¹BA⁻¹, C⁻¹, −C⁻¹f/e, 1/e)`, `A⁻¹` and `C⁻¹` by `Mat3::inverse_adj`.
    ///
    /// # Domain
    ///
    /// `e ≠ 0`, checked by `debug_assert!` ("no lane is zero", so a NaN Jacobian returns NaN):
    /// `ad`'s `e = 0` is the structural singularity. `det A`, `det C` nonzero and normal, as
    /// `Mat3::inverse_adj` states, unasserted: `C = s R` or `V` has `det C` of order `e^{3σ}`, which
    /// under- or overflows only outside `NUMERICS.md` §12's `|σ| < ln MAX`. A release build divides
    /// (D11).
    #[inline]
    fn inverse(&self) -> Self {
        let (ai, _) = self.a.inverse_adj();
        let (ci, _) = self.c.inverse_adj();
        debug_assert!(!self.e.abs().le(S::zero()).any(), "Sim3Jac::inverse: e = 0");
        let ei = S::one() / self.e;
        Self {
            a: ai,
            b: -(ci * self.b * ai),
            c: ci,
            f: -(ci * self.f).scale(ei),
            e: ei,
        }
    }
    #[inline]
    fn neg(&self) -> Self {
        Self {
            a: -self.a,
            b: -self.b,
            c: -self.c,
            f: -self.f,
            e: -self.e,
        }
    }
    #[inline]
    fn apply(&self, t: &Sim3Tangent<S>) -> Sim3Tangent<S> {
        Sim3Tangent {
            phi: self.a * t.phi,
            rho: (self.b * t.phi + self.c * t.rho) + self.f.scale(t.sigma),
            sigma: self.e * t.sigma,
        }
    }
    /// `[Aᵀφ + Bᵀρ; Cᵀρ; f·ρ + eσ]`.
    #[inline]
    fn apply_transpose(&self, t: &Sim3Tangent<S>) -> Sim3Tangent<S> {
        Sim3Tangent {
            phi: self.a.transpose() * t.phi + self.b.transpose() * t.rho,
            rho: self.c.transpose() * t.rho,
            sigma: self.f.dot(t.rho) + self.e * t.sigma,
        }
    }
    #[inline]
    fn write_dense(&self, out: &mut StridedMut<'_, S>) {
        debug_assert!(
            out.rows() == 7 && out.cols() == 7,
            "Jac::write_dense: the view is not DOF x DOF"
        );
        for c in 0..7 {
            for r in 0..7 {
                let v = match (r, c) {
                    (0..=2, 0..=2) => self.a.get(r, c),
                    (3..=5, 0..=2) => self.b.get(r - 3, c),
                    (3..=5, 3..=5) => self.c.get(r - 3, c - 3),
                    (3..=5, 6) => self.f.0[r - 3],
                    (6, 6) => self.e,
                    _ => S::zero(),
                };
                out.set(r, c, v);
            }
        }
    }
    /// `J Σ Jᵀ`: `M = J Σ` by columns, then row `r` of `M Jᵀ` is `J` applied to row `r` of `M`
    /// (`SEn3Jac::sandwich`'s order), each through the structure, skipping its zeros.
    #[inline]
    fn sandwich<const D: usize>(&self, cov: &Matrix<S, D, D>) -> Matrix<S, D, D> {
        const { assert!(D == <Sim3Tangent<S> as Tangent<S>>::DOF) };
        let m = Matrix::from_cols(array::from_fn(|c| apply7::<S, D>(self, &cov.col(c))));
        Matrix::from_rows(array::from_fn(|r| apply7::<S, D>(self, &m.row(r))))
    }
}
