//! SO(3): the rotation group on the unit quaternion (`docs/PHASE3.md` §4).
//!
//! `Exp` and `Log` are `NUMERICS.md` §3.1 and §3.2 over the grouped kernels of
//! [`crate::coeffs`], which is what makes this the kernel's first consumer: `exp_coeffs` answers
//! `(k, cos θ/2)` in **one** `branch`, as §4's group rule asks, and `log_ratio` takes the swept
//! `r` switch on `s = n²/w²`. The conformance subject that scores those two kernels over
//! `so3_exp` and `so3_log` is `xtask`'s; this type is the shipped path it scores.
//!
//! `Jac = Mat3<S>`: SO(3)'s adjoint *is* a rotation matrix (`NUMERICS.md` §3.5,
//! `Ad_R = R`), so there is no structure to exploit below a dense `3 x 3` and `0005`'s structured
//! type is the matrix itself.

use crate::coeffs::{exp_coeffs, jr_coeffs, jr_inv_coeff, log_ratio};
use crate::quat::Quat;
use crate::side::Side;
use crate::traits::{tie_dof, Jac, LieGroup, Tangent};
use core::array;
use core::ops::Mul;
use helicoid_linalg::{hat, Blend, Mask, Mat3, Matrix, Real, StridedMut, Vec3, Vector};

/// A rotation, stored as a unit quaternion (Hamilton, `w` first, active — `0002`).
///
/// The group law is `Mul`, never `Add`:
///
/// ```
/// use helicoid::{LieGroup, SO3, SO3Tangent};
/// use helicoid_linalg::Vector;
/// let r = SO3::<f64>::exp(&SO3Tangent { phi: Vector([0.0, 0.0, core::f64::consts::FRAC_PI_2]) });
/// let p = r.act(Vector([1.0, 0.0, 0.0]));
/// assert!((p.0[1] - 1.0).abs() < 1e-15);
/// ```
#[repr(transparent)]
#[derive(Clone, Copy, Debug)]
pub struct SO3<S>(Quat<S>);

/// A tangent vector of [`SO3`]: the rotation vector `phi`, in dense order (`0002`).
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct SO3Tangent<S> {
    /// The rotation vector `φ = θ n̂`; `write_dense` writes its components in index order.
    pub phi: Vec3<S>,
}

/// `x·x + y·y + z·z` left to right: `θ²` and `n²` are dot products (`NUMERICS.md` §2.1), and this
/// order is the one the kernels and the corpus were scored against.
#[inline]
fn norm_sq<S: Real>(v: Vec3<S>) -> S {
    let [x, y, z] = v.0;
    (x * x + y * y) + z * z
}

/// `m` applied to a `D`-vector, for the `sandwich` whose `D` is `3` by its own assertion. Written
/// out because a `Matrix<S, D, D>` column is not a `Vec3<S>` to the type checker, as `SEn3Jac`
/// does the same thing for the same reason; the sum is left to right, as `Mul<Vector>` forms it,
/// so every entry is bit-identical to the dense product.
#[inline]
fn apply_rows<S: Real, const D: usize>(m: &Mat3<S>, v: &Vector<S, D>) -> Vector<S, D> {
    let rows = [m.row(0).0, m.row(1).0, m.row(2).0];
    Vector(array::from_fn(|r| {
        let a = rows[r];
        (a[0] * v.0[0] + a[1] * v.0[1]) + a[2] * v.0[2]
    }))
}

impl<S: Real> SO3<S> {
    /// The rotation of a quaternion the caller guarantees is unit, within `NUMERICS.md` §3.6's
    /// bound; `Quat::from_wxyz_unchecked` carries the `debug_assert!`.
    #[inline]
    pub fn from_quat_unchecked(q: Quat<S>) -> Self {
        Self(q)
    }

    /// The rotation of any non-zero quaternion, normalized (`0027` decision 4).
    #[inline]
    pub fn from_quat_normalized(q: Quat<S>) -> Self {
        Self(Quat::from_wxyz_normalized(q.w, q.x, q.y, q.z))
    }

    /// The stored quaternion.
    #[inline]
    pub fn quat(&self) -> Quat<S> {
        self.0
    }

    /// `R(q)`, the rotation matrix (`NUMERICS.md` §1).
    #[inline]
    pub fn to_matrix(&self) -> Mat3<S> {
        self.0.to_matrix()
    }

    /// The rotation of a `3 x 3` matrix by Shepperd's method (`NUMERICS.md` §3.4): the largest of
    /// `{tr R, R00, R11, R22}` through nested `branch`es, then normalize. **Closed form, never
    /// iterative** — locus-tag's `quat_from_so3` regression is an iteration that does not
    /// terminate on degenerate input.
    ///
    /// Each arm takes its own square root at a *safe* argument (`0003` item 3): the four
    /// candidates sum to `4`, so the selected one is at least `1`, and an unselected arm is handed
    /// `1` instead of a possibly negative value. So no lane ever evaluates `sqrt` of a negative
    /// number or divides by zero, whatever the input matrix.
    ///
    /// # Domain
    ///
    /// `r` is a rotation matrix. For anything else this returns the nearest quaternion its pivot
    /// admits, which is not a projection onto SO(3); §3.4 says project with `svd3` first. It does
    /// not panic, in debug or release: an entry large enough to overflow a pivot candidate gives a
    /// non-finite or NaN quaternion, as `Quat::from_wxyz_unchecked` does -- garbage in, garbage
    /// out. That is why the normalization below is written out rather than taken from
    /// `Quat::from_wxyz_normalized`, whose `debug_assert!` would turn such a matrix, and an
    /// all-NaN one, into a panic. `SO3::exp` and `coeffs`'s `nonnegative` read NaN the same way.
    #[inline]
    pub fn from_matrix(r: &Mat3<S>) -> Self {
        // `get` per entry, each read once: `row(i)` builds and discards a whole `Vec3` per scalar,
        // and the diagonal is wanted four times over (the trace and the three axis candidates).
        let (d0, d1, d2) = (r.get(0, 0), r.get(1, 1), r.get(2, 2));
        let (two, half) = (S::lit(2.0), S::lit(0.5));
        let trace = (d0 + d1) + d2;
        let c = [
            S::one() + trace,
            (S::one() + d0 - d1) - d2,
            (S::one() - d0 + d1) - d2,
            (S::one() - d0 - d1) + d2,
        ];
        // `lt` only: `b.lt(a)` is `a > b`, and `Real` has no `PartialOrd` (`0003`).
        let best01 = c[1].lt(c[0]);
        let best23 = c[3].lt(c[2]);
        let first = S::select(best23, c[2], c[3]).lt(S::select(best01, c[0], c[1]));
        let pick = [
            first.and(best01),
            first.and(best01.not()),
            first.not().and(best23),
            first.not().and(best23.not()),
        ];
        // The three off-diagonal differences and sums every arm reads.
        let (m01, m10) = (r.get(0, 1), r.get(1, 0));
        let (m02, m20) = (r.get(0, 2), r.get(2, 0));
        let (m12, m21) = (r.get(1, 2), r.get(2, 1));
        let (dx, dy, dz) = (m21 - m12, m02 - m20, m10 - m01);
        let (sx, sy, sz) = (m01 + m10, m02 + m20, m12 + m21);
        // `d >= 1` when selected and exactly `1` when not, so `two * d` never divides by zero.
        let arm = |i: usize| {
            let d = S::select(pick[i], c[i], S::one()).sqrt();
            (half * d, two * d)
        };
        let q = S::branch(
            first,
            || {
                S::branch(
                    best01,
                    || {
                        let (p, t) = arm(0);
                        [p, dx / t, dy / t, dz / t]
                    },
                    || {
                        let (p, t) = arm(1);
                        [dx / t, p, sx / t, sy / t]
                    },
                )
            },
            || {
                S::branch(
                    best23,
                    || {
                        let (p, t) = arm(2);
                        [dy / t, sx / t, p, sz / t]
                    },
                    || {
                        let (p, t) = arm(3);
                        [dz / t, sy / t, sz / t, p]
                    },
                )
            },
        );
        // The same four divisions `Quat::from_wxyz_normalized` performs, in the same order and on
        // the same `norm_sq`, so an in-domain matrix gives the identical quaternion -- without its
        // `debug_assert!`, for the reason the *Domain* note gives.
        let n = Quat {
            w: q[0],
            x: q[1],
            y: q[2],
            z: q[3],
        }
        .norm_sq()
        .sqrt();
        Self(Quat {
            w: q[0] / n,
            x: q[1] / n,
            y: q[2] / n,
            z: q[3] / n,
        })
    }

    /// `R v` by the quaternion sandwich (`NUMERICS.md` §3.3):
    /// `v' = v + 2w (u × v) + 2 u × (u × v)`, in that order.
    ///
    /// # Domain
    ///
    /// A unit `q`. The formula equals `R(q)v` only there: it expands to
    /// `(1 − 2‖u‖²)v + 2w(u×v) + 2u(u·v)`, where `R(q)` has `w² − ‖u‖²` in place of `1 − 2‖u‖²`.
    /// Off the unit sphere the three candidate readings are measurably different, and the corpus's
    /// `so3_act` reference takes the **normalized** one: at the `q:nonunit` stratum
    /// (`η = ‖q‖² − 1 = 2^-45`) this formula scores `510.3 u` against it and `to_matrix() * v`,
    /// the scaled rotation of §1, scores `257.8 u` — which are `2·2^-45/u` and `2^-45/u` exactly.
    /// Which reading `so3_act` names is an open item of `docs/maths/index.md`; this method is §3.3
    /// as written and normalizes nothing.
    #[inline]
    pub fn act(&self, v: Vec3<S>) -> Vec3<S> {
        let two = S::lit(2.0);
        let u = Vector([self.0.x, self.0.y, self.0.z]);
        let uv = u.cross(v);
        (v + uv.scale(two * self.0.w)) + u.cross(uv).scale(two)
    }

    /// `R` applied to every point in place, forming `R(q)` once (`NUMERICS.md` §3.3).
    ///
    /// Its reference twin is the per-point [`act`](SO3::act) (`NUMERICS.md` §14, tolerance `3 u`),
    /// compared against it by `act_many_matches_reference` as D6 requires. The two are *not*
    /// bit-identical: `R(q)` rounds the nine matrix entries before any point is touched, where
    /// `act` rounds the sandwich per point.
    #[inline]
    pub fn act_many(&self, pts: &mut [Vec3<S>]) {
        let r = self.to_matrix();
        for p in pts.iter_mut() {
            *p = r * *p;
        }
    }

    /// One Newton step back onto the unit sphere (`NUMERICS.md` §3.6).
    #[inline]
    pub fn renormalize(&mut self) {
        self.0.renormalize();
    }
}

impl<S: Real> Mul for SO3<S> {
    type Output = Self;
    /// The Hamilton product: `a * b` is `R_a_x · R_x_b` (`0002`).
    #[inline]
    fn mul(self, o: Self) -> Self {
        Self(self.0 * o.0)
    }
}

impl<S: Real> Blend<S> for SO3<S> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self(Quat::blend(m, t.0, f.0))
    }
}

impl<S: Real> Blend<S> for SO3Tangent<S> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self {
            phi: Vector::blend(m, t.phi, f.phi),
        }
    }
}

impl<S: Real> Tangent<S> for SO3Tangent<S> {
    const DOF: usize = 3;
    #[inline]
    fn zero() -> Self {
        Self {
            phi: Vector([S::zero(); 3]),
        }
    }
    #[inline]
    fn add(&self, o: &Self) -> Self {
        Self {
            phi: self.phi + o.phi,
        }
    }
    #[inline]
    fn sub(&self, o: &Self) -> Self {
        Self {
            phi: self.phi - o.phi,
        }
    }
    #[inline]
    fn neg(&self) -> Self {
        Self { phi: -self.phi }
    }
    #[inline]
    fn scale(&self, k: S) -> Self {
        Self {
            phi: self.phi.scale(k),
        }
    }
    #[inline]
    fn dot_acc(&self, o: &Self, acc: S) -> S {
        // `Vector::dot` would start its own sum, which is the same answer only for `acc = +0`.
        self.phi
            .0
            .iter()
            .zip(o.phi.0)
            .fold(acc, |a, (&x, y)| a + x * y)
    }
    #[inline]
    fn write_dense(&self, out: &mut [S]) {
        debug_assert!(out.len() == 3, "Tangent::write_dense: wrong length");
        match out.first_chunk_mut::<3>() {
            Some(o) => *o = self.phi.0,
            None => {
                for (o, v) in out.iter_mut().zip(self.phi.0) {
                    *o = v;
                }
            }
        }
    }
    #[inline]
    fn read_dense(src: &[S]) -> Self {
        debug_assert!(src.len() == 3, "Tangent::read_dense: wrong length");
        // A short `src` reads NaN, not `+0`: see `RnTangent::read_dense` for why (D11).
        let phi = match src.first_chunk::<3>() {
            Some(s) => Vector(*s),
            None => Vector(array::from_fn(|i| {
                src.get(i).copied().unwrap_or_else(|| S::zero() / S::zero())
            })),
        };
        Self { phi }
    }
}

impl<S: Real> Jac<S, SO3Tangent<S>> for Mat3<S> {
    #[inline]
    fn identity() -> Self {
        Matrix::identity()
    }
    #[inline]
    fn mul(&self, o: &Self) -> Self {
        *self * *o
    }
    /// The algebraic inverse, by `Mat3::inverse_adj`: the adjugate over the determinant (`0005` —
    /// not a numerical inverse, and not a transpose, because `ad` and `J_r` are not orthogonal).
    ///
    /// `inverse_adj` divides by `det` itself and returns `(J⁻¹, det)`, so nothing is scaled here;
    /// scaling the first element again is a `1/det²` that `laws::jac_dense_order` catches at
    /// `1.5e14 u`.
    ///
    /// # Domain
    ///
    /// `det != 0`, checked by `debug_assert!`; `ad_φ = W` is singular for every `φ` and has no
    /// inverse. A release build divides and returns `±inf` or NaN (D11).
    ///
    /// The check reads "no lane has `det = 0`", not "every lane has `|det| > 0`", for the reason
    /// `coeffs`'s `nonnegative` does: the two differ only on NaN, and a NaN determinant is a NaN
    /// matrix, whose inverse is NaN — the answer, not a panic.
    #[inline]
    fn inverse(&self) -> Self {
        let (inv, det) = self.inverse_adj();
        debug_assert!(
            !det.abs().le(S::zero()).any(),
            "Jac::inverse: singular 3 x 3"
        );
        inv
    }
    #[inline]
    fn neg(&self) -> Self {
        -*self
    }
    #[inline]
    fn apply(&self, t: &SO3Tangent<S>) -> SO3Tangent<S> {
        SO3Tangent { phi: *self * t.phi }
    }
    #[inline]
    fn apply_transpose(&self, t: &SO3Tangent<S>) -> SO3Tangent<S> {
        SO3Tangent {
            phi: self.transpose() * t.phi,
        }
    }
    #[inline]
    fn write_dense(&self, out: &mut StridedMut<'_, S>) {
        debug_assert!(
            out.rows() == 3 && out.cols() == 3,
            "Jac::write_dense: the view is not DOF x DOF"
        );
        for c in 0..3 {
            for (r, v) in self.col(c).0.iter().enumerate() {
                out.set(r, c, *v);
            }
        }
    }
    #[inline]
    fn sandwich<const D: usize>(&self, cov: &Matrix<S, D, D>) -> Matrix<S, D, D> {
        const { assert!(D == <SO3Tangent<S> as Tangent<S>>::DOF) };
        // `M = J Σ` by columns, then `M Jᵀ` by rows: row `r` of `M Jᵀ` is `J` applied to row `r`
        // of `M`, so no transpose and no `D x D` scratch is formed (`SEn3Jac::sandwich`).
        let m = Matrix::from_cols(array::from_fn(|c| apply_rows::<S, D>(self, &cov.col(c))));
        Matrix::from_rows(array::from_fn(|r| apply_rows::<S, D>(self, &m.row(r))))
    }
}

impl<S: Real> LieGroup<S> for SO3<S> {
    type Tangent = SO3Tangent<S>;
    type Jac = Mat3<S>;
    const DOF: usize = 3;

    #[inline]
    fn identity() -> Self {
        const { tie_dof::<S, Self>() };
        Self(Quat::from_wxyz_unchecked(
            S::one(),
            S::zero(),
            S::zero(),
            S::zero(),
        ))
    }
    /// The conjugate: `R⁻¹ = Rᵀ` and `q⁻¹ = q*` for a unit `q`.
    #[inline]
    fn inverse(&self) -> Self {
        Self(self.0.conjugate())
    }
    /// `q = (cos(θ/2), k(θ) φ)` (`NUMERICS.md` §3.1), both coefficients from one grouped
    /// `branch` of `exp_coeffs`.
    #[inline]
    fn exp(tau: &SO3Tangent<S>) -> Self {
        let (k, cos_half) = exp_coeffs(norm_sq(tau.phi));
        let [x, y, z] = tau.phi.0;
        // Built field by field, not through `Quat::from_wxyz_unchecked`. That constructor's
        // `debug_assert!` is for a *caller* handing in a quaternion it claims is unit; `Exp`'s
        // output is unit by construction, to the rounding its own strata record, so the check has
        // nothing to add here. It would also turn a NaN tangent into a panic, where a value
        // function owes its caller NaN — which is what `laws::dual_value_is_plain_value` asks for.
        Self(Quat {
            w: cos_half,
            x: k * x,
            y: k * y,
            z: k * z,
        })
    }
    /// `NUMERICS.md` §3.2: the `copysign` flip, then `φ = r(n², w) u`.
    ///
    /// Exactly scale-invariant in `q`, as `atan2` and `u/n` are, so a `q` that is unit only within
    /// rounding returns the `Log` of its normalization. At `w = +0` nothing flips, so `Log` is a
    /// function of the **quaternion**, not of the rotation: `q` and `-q` return `±π n̂`.
    #[inline]
    fn log(&self) -> SO3Tangent<S> {
        let flip = S::one().copysign(self.0.w);
        let (w, x, y, z) = (
            flip * self.0.w,
            flip * self.0.x,
            flip * self.0.y,
            flip * self.0.z,
        );
        let r = log_ratio((x * x + y * y) + z * z, w);
        SO3Tangent {
            phi: Vector([r * x, r * y, r * z]),
        }
    }
    /// `Ad_R = R` (`NUMERICS.md` §3.5).
    #[inline]
    fn adjoint(&self) -> Mat3<S> {
        self.to_matrix()
    }
    /// `ad_φ = W = [φ]_×` (`NUMERICS.md` §3.5).
    #[inline]
    fn ad(tau: &SO3Tangent<S>) -> Mat3<S> {
        hat(tau.phi)
    }
    /// `J_r = I − aW + bW²` (`NUMERICS.md` §3.5), `(a, b)` from one grouped `branch`.
    ///
    /// `jl` is **not** overridden: `J_l(φ) = J_r(−φ)` is the trait's provided body, and it is
    /// `I + aW + bW²` to the bit — `a` and `b` are functions of `θ²`, and `(−W)² = W²` entry by
    /// entry in the same products and the same order. One code path, and §14's twin stays the
    /// provided one.
    #[inline]
    fn jr(tau: &SO3Tangent<S>) -> Mat3<S> {
        let w = hat(tau.phi);
        let (a, b) = jr_coeffs(norm_sq(tau.phi));
        // `W²` through the generic product, not the closed `φφᵗ − θ²I`. The closed form costs 6
        // multiplies against 27 and lowers `so3_jr`'s worst row from 4.097 to 3.439 `u`, but it is
        // a different rounding (13.8% of entries differ over 20 000 samples) and the corpus says
        // it is worse where it is not better: 6 of 28 `so3_jr` strata regress, up to 1.34x, and
        // `so3_jr_inv` keeps its 2.112 maximum while **12 of 28** strata regress, up to 1.39x.
        // Nothing has asked for the arithmetic yet — `PHASE3.md` §11's benches are owed — so the
        // trade is not taken, and `0006` says the bar is the max, per stratum, never a mean.
        (Matrix::identity() + w.scale(-a)) + (w * w).scale(b)
    }
    /// `J_r⁻¹ = I + W/2 + cW²` (`NUMERICS.md` §3.5).
    ///
    /// # Domain
    ///
    /// `θ < 2π` (`NUMERICS.md` §12): `J_r` is singular at `2π` and `c` has a pole there. Nothing
    /// is asserted — `θ` is not a stored field and the check would cost a norm on every call.
    #[inline]
    fn jr_inv(tau: &SO3Tangent<S>) -> Mat3<S> {
        let w = hat(tau.phi);
        let c = jr_inv_coeff(norm_sq(tau.phi));
        (Matrix::identity() + w.scale(S::lit(0.5))) + (w * w).scale(c)
    }
    /// `(Ad_Exp(τ)⁻¹, J_r(τ))` (`NUMERICS.md` §2.3). `Ad` is a rotation here, so its inverse is
    /// the transpose, not `Jac::inverse`'s adjugate.
    #[inline]
    fn rplus_jacobians(&self, tau: &SO3Tangent<S>) -> (Mat3<S>, Mat3<S>) {
        (Self::exp(tau).to_matrix().transpose(), Self::jr(tau))
    }
    /// `(Ad_Exp(τ), J_l(τ))` (`NUMERICS.md` §2.3).
    #[inline]
    fn lplus_jacobians(&self, tau: &SO3Tangent<S>) -> (Mat3<S>, Mat3<S>) {
        (Self::exp(tau).to_matrix(), Self::jl(tau))
    }
    /// `(J_r⁻¹(τ), −J_l⁻¹(τ))` at `τ = self ⊖_R base` (`NUMERICS.md` §2.3), from **one**
    /// `jr_inv` and a transpose.
    ///
    /// `J_l(φ) = J_r(φ)ᵗ` and `J_l⁻¹(φ) = J_r⁻¹(φ)ᵗ` **bit for bit**, not just mathematically:
    /// `J_r = I − aW + bW²` with `Wᵗ = −W` gives `J_rᵗ = I + aW + bW²`, and in floating point
    /// `hat(−φ)` is exactly `hat(φ)ᵗ`, `a` and `b` are functions of `θ²` alone, and `W²` as
    /// `Matrix::mul` forms it is exactly symmetric (its `(i, j)` and `(j, i)` sums are the same
    /// products in the same index order). Pinned by `jl_is_jr_transposed_to_the_bit`.
    ///
    /// This row and `lminus_jacobians` are the crate's hottest Jacobian path — one per residual
    /// per solver iteration — and computing both inverses independently paid a second `norm_sq`,
    /// a second grouped coefficient `branch` (two `sqrt` and two `sin_cos` on the exact arm), a
    /// second `hat` and a second 27-multiply `3 x 3` product for a matrix a transpose already has.
    #[inline]
    fn rminus_jacobians(&self, base: &Self) -> (Mat3<S>, Mat3<S>) {
        let jri = Self::jr_inv(&self.rminus(base));
        (jri, -jri.transpose())
    }
    /// `(J_l⁻¹(τ), −J_r⁻¹(τ))` at `τ = self ⊖_L base` (`NUMERICS.md` §2.3), from one `jr_inv`
    /// and a transpose; see [`rminus_jacobians`](SO3::rminus_jacobians) for why that is exact.
    #[inline]
    fn lminus_jacobians(&self, base: &Self) -> (Mat3<S>, Mat3<S>) {
        let jri = Self::jr_inv(&self.lminus(base));
        (jri.transpose(), -jri)
    }
    /// Right `(Ad_Y⁻¹, I)`, left `(I, Ad_X)` (`NUMERICS.md` §2.3), selected by `Sd::IS_RIGHT` at
    /// compile time.
    #[inline]
    fn compose_jacobians<Sd: Side>(&self, rhs: &Self) -> (Mat3<S>, Mat3<S>) {
        match Sd::IS_RIGHT {
            true => (rhs.to_matrix().transpose(), Matrix::identity()),
            false => (Matrix::identity(), self.to_matrix()),
        }
    }
    /// Right `−Ad_X`, left `−Ad_X⁻¹` (`NUMERICS.md` §2.3).
    #[inline]
    fn inverse_jacobian<Sd: Side>(&self) -> Mat3<S> {
        match Sd::IS_RIGHT {
            true => -self.to_matrix(),
            false => -self.to_matrix().transpose(),
        }
    }
}
