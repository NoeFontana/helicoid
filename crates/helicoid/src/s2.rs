//! S² and its frozen Householder chart (`docs/PHASE5.md` §2, `NUMERICS.md` §8, `0063`, `0067`).
//!
//! A point is a unit `n ∈ ℝ³`; the chart at `n` is the basis `B = [b₁ b₂]` of `n^⊥` computed
//! once in [`Chart::at`](crate::Chart::at) and reused by every call: the basis is discontinuous
//! at `n_z = 0`, as every global basis is somewhere (`docs/maths/charts.md` CH.11), so a chart
//! that recomputed it from a moving point would change what its coordinates mean (`0012`).

use crate::chart::{Chart, Manifold};
use crate::coeffs::{exp_coeffs, jr_coeffs, log_ratio};
use crate::quat::unit_tolerance;
use crate::rn::RnTangent;
use crate::traits::{Jac, Tangent};
use core::array;
use helicoid_linalg::{Blend, Mask, Mat2, Matrix, Real, StridedMut, Vec3, Vector};

/// A point of the unit sphere S², held unit (`0063` decision 2).
///
/// Built by one of three entry points, as the quaternion is:
/// [`from_vec_unchecked`](S2::from_vec_unchecked) vouches,
/// [`from_vec_normalized`](S2::from_vec_normalized) divides by the norm, and
/// [`renormalize`](S2::renormalize) repairs drift with one Newton step.
#[repr(transparent)]
#[derive(Clone, Copy, Debug)]
pub struct S2<S>(pub(crate) Vec3<S>);

impl<S: Real> S2<S> {
    /// The point `v`, taken as it is.
    ///
    /// # Domain
    ///
    /// `| ‖v‖² - 1 | <= 2^-40` (`f64`, `Dual<f64, N>`) or `2^-16` (`f32`) (`NUMERICS.md` §12),
    /// checked by `debug_assert!` (a NaN fails it). Release builds check nothing: garbage in,
    /// garbage out. One [`renormalize`](S2::renormalize) step repairs a point inside it.
    #[inline]
    pub fn from_vec_unchecked(v: Vec3<S>) -> Self {
        debug_assert!(
            (v.norm_sq() - S::one())
                .abs()
                .le(unit_tolerance::<S>())
                .all(),
            "S2::from_vec_unchecked: | |v|^2 - 1 | exceeds 2^-40 (f64) or 2^-16 (f32)"
        );
        Self(v)
    }

    /// `v / ‖v‖`: a normalization, for a direction whose norm the caller cannot vouch for
    /// (`0027`). One `sqrt` and one division per component.
    ///
    /// # Domain
    ///
    /// `‖v‖²` normal: nonzero, and the components within about `10^±154` (`f64`) or `10^±19`
    /// (`f32`), `debug_assert!`ed on the result. Release builds check nothing: at `v = 0` every
    /// component is `0/0`, and above the overflow the result is zero.
    #[inline]
    pub fn from_vec_normalized(v: Vec3<S>) -> Self {
        let n = v.norm();
        let p = Vector(v.0.map(|x| x / n));
        debug_assert!(
            (p.norm_sq() - S::one())
                .abs()
                .le(unit_tolerance::<S>())
                .all(),
            "S2::from_vec_normalized: |v|^2 is zero, not finite, or NaN"
        );
        Self(p)
    }

    /// One Newton step, `n ← n (3 - ‖n‖²)/2`, in place: the quaternion's step on three
    /// components (`NUMERICS.md` §3.6, `API.md` R2).
    ///
    /// A normalization for `| ‖n‖² - 1 | <= 2^-26.29` (`f64`) or `2^-11.79` (`f32`): one step lands
    /// within `5u` of the sphere there (measured, `0063`). Beyond it, only a better guess.
    ///
    /// # Domain
    ///
    /// None is asserted; the step is defined for every input, as the quaternion's is.
    #[inline]
    pub fn renormalize(&mut self) {
        let k = (S::lit(3.0) - self.0.norm_sq()) * S::lit(0.5);
        self.0 = self.0.scale(k);
    }

    /// The unit vector.
    #[inline]
    pub fn vec(&self) -> Vec3<S> {
        self.0
    }
}

impl<S: Real> Blend<S> for S2<S> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self(Vector::blend(m, t.0, f.0))
    }
}

/// The chart of [`S2`] frozen at a base `n` (`NUMERICS.md` §8): coordinates `δ ∈ ℝ²` in the
/// Householder basis `B = [b₁ b₂]` of `n^⊥`.
///
/// `retract(δ) = Exp(Bδ) n`; `local(m) = Bᵀ Log(E)`, `E` the minimal rotation taking `n` to `m`.
/// `δ` is a **rotation vector** `Bδ`, a quarter turn from the displacement it causes
/// (`docs/maths/charts.md` CH.9(a)). The Jacobians are `Mat2` and read the perturbed point in its
/// own frozen frame (CH.10), so they carry the frame change `ς ς_m` between the two bases.
#[derive(Clone, Copy, Debug)]
pub struct S2Chart<S> {
    base: S2<S>,
    pub(crate) b1: Vec3<S>,
    pub(crate) b2: Vec3<S>,
}

/// `ς`: `+1` unless `n_z < 0`, a comparison, so `±0` and NaN give `+1` (`0063` decision 1).
#[inline]
fn sign<S: Real>(n: Vec3<S>) -> S {
    S::select(n.0[2].lt(S::zero()), -S::one(), S::one())
}

/// `(H e_x, H e_y)`, `H = I - 2ννᵀ/(νᵀν)`, `ν = n + ς e_z` (CH.8). The `H` form, with `νᵀν`
/// computed rather than taken as `2(1 + |n_z|)`: it keeps `BᵀB = I` to rounding for a point off
/// the sphere, where the closed form CH.8(d) is off by twice the drift in `BᵀB` and `Bᵀn` (CH.8(e)).
#[inline]
fn basis<S: Real>(n: Vec3<S>) -> (Vec3<S>, Vec3<S>) {
    let [x, y, z] = n.0;
    let nu = Vector([x, y, z + sign(n)]);
    let k = S::lit(2.0) / nu.norm_sq();
    let (kx, ky) = (k * x, k * y);
    (
        Vector([S::one() - kx * x, -(kx * y), -(kx * nu.0[2])]),
        Vector([-(ky * x), S::one() - ky * y, -(ky * nu.0[2])]),
    )
}

impl<S: Real> S2Chart<S> {
    /// `(retract(δ), θ²)`: `cos θ n + ς (sin θ/θ)(δ₂ b₁ - δ₁ b₂)`, which is `Exp(Bδ) n` since
    /// `Bδ × n = -ς B K δ` (CH.9(a)). `sin θ/θ = 2k cos(θ/2)` and `cos θ = 1 - 2k²θ²` from the one
    /// `exp_coeffs`, finite at `δ = 0`. [`retract_jacobian`](Chart::retract_jacobian) reads the
    /// point from here, so the frame of `m` it is expressed in is the one `retract` returns.
    #[inline]
    fn step(&self, delta: &RnTangent<S, 2>) -> (Vec3<S>, S) {
        let [d1, d2] = delta.rho.0;
        let z = delta.rho.norm_sq();
        let (k, cos_half) = exp_coeffs(z);
        let cos = S::one() - S::lit(2.0) * (k * k * z);
        let g = sign(self.base.0) * (S::lit(2.0) * k * cos_half);
        let n = self.base.0;
        let m = Vector(array::from_fn(|i| {
            cos * n.0[i] + g * (d2 * self.b1.0[i] - d1 * self.b2.0[i])
        }));
        (m, z)
    }

    /// `local(other)` and what its Jacobian reuses (CH.9(b)): `δ = Bᵀp`, `p = (α/s)(n × m)` the
    /// rotation vector of the minimal rotation (`‖p‖ = α`), `w = n·m` and the ratio `α/s`.
    ///
    /// `α/s = ½ r(s², w)` from `log_ratio`, whose mask on `s²/w²` takes the series arm only for
    /// `w > 0` (CH.9(c), `0063` decision 3). `Bᵀ(α/s · n × m)` and not `ς (α/s) K Bᵀ m`: the first
    /// measured 1.4–1.9 `u/(π - α)` near the antipode, the second 2.3–3.8 (CH.11(e)). At `m = -n`
    /// (`s² = 0`, `w < 0`) the exact arm would divide `2π` by `0`; the safe argument `s² = 1` makes
    /// the answer `0` there, the finite value `NUMERICS.md` §12 promises outside the domain.
    #[inline]
    fn log_parts(&self, other: &S2<S>) -> LogParts<S> {
        let (n, m) = (self.base.0, other.0);
        let v = n.cross(m);
        let w = n.dot(m);
        let s2 = v.norm_sq();
        let antipode = s2.le(S::zero()).and(w.le(S::zero()));
        let ratio = S::lit(0.5) * log_ratio(S::select(antipode, S::one(), s2), w);
        let p = v.scale(ratio);
        let delta = RnTangent {
            rho: Vector([self.b1.dot(p), self.b2.dot(p)]),
        };
        LogParts { delta, w, ratio, p }
    }

    /// `local_jacobian` from [`log_parts`](Self::log_parts): `A(δ)⁻¹ Z` (CH.10(b)),
    /// `A⁻¹ = (α/s)(I - bδδᵀ)`, `Z_ij = (E b_i)·b'_j` with `[b'₁ b'₂]` the basis at `m`.
    ///
    /// `E b_i = w b_i - (b_i·m) n + a δ_i p` is Rodrigues for `E = Exp(p)` reduced by `b_i ⊥ n`
    /// (`(n × m) × b_i = -(b_i·m) n`, `p·b_i = δ_i`, `1 - a‖p‖² = w`). It never divides by
    /// `1 + w ≈ (π - α)²/2`, as CH.10(c)'s `(v, w)` form does, which measured `1.5e9 u` at
    /// `π - α = 10⁻⁴` against `9.6 u` for `Exp` (CH.10, *Checked*). `a`, `b` are `jr_coeffs` at
    /// `‖δ‖² = α²`.
    #[inline]
    fn local_jac(&self, other: &S2<S>, l: &LogParts<S>) -> Mat2<S> {
        let (n, m) = (self.base.0, other.0);
        let [d1, d2] = l.delta.rho.0;
        let (a, b) = jr_coeffs(l.delta.rho.norm_sq());
        let (c1, c2) = basis(m);
        let e = |bi: Vec3<S>, di: S| {
            let (h, ad) = (bi.dot(m), a * di);
            Vector(array::from_fn(|k| {
                (l.w * bi.0[k] - h * n.0[k]) + ad * l.p.0[k]
            }))
        };
        let (e1, e2) = (e(self.b1, d1), e(self.b2, d2));
        // `Z` by rows: row `i` is `E b_i` read in the frame at `m`.
        let zr = [[e1.dot(c1), e1.dot(c2)], [e2.dot(c1), e2.dot(c2)]];
        let u = [d1 * zr[0][0] + d2 * zr[1][0], d1 * zr[0][1] + d2 * zr[1][1]];
        let bd = [b * d1, b * d2];
        Matrix::from_rows(array::from_fn(|r| {
            Vector(array::from_fn(|c| l.ratio * (zr[r][c] - bd[r] * u[c])))
        }))
    }
}

/// What [`S2Chart::log_parts`] computes once for `local` and `local_jacobian`.
#[derive(Clone, Copy)]
struct LogParts<S> {
    delta: RnTangent<S, 2>,
    w: S,
    ratio: S,
    p: Vec3<S>,
}

impl<S: Real> Chart<S, S2<S>> for S2Chart<S> {
    type Tangent = RnTangent<S, 2>;
    type Jac = Mat2<S>;
    /// The basis at `base`, its one computation (`NUMERICS.md` §8).
    #[inline]
    fn at(base: &S2<S>) -> Self {
        let (b1, b2) = basis(base.0);
        Self {
            base: *base,
            b1,
            b2,
        }
    }
    #[inline]
    fn base(&self) -> S2<S> {
        self.base
    }
    /// `Exp(Bδ) n`, for every `δ` (CH.9(a)). The result is as unit as the base, to a few `u`.
    #[inline]
    fn retract(&self, delta: &RnTangent<S, 2>) -> S2<S> {
        S2(self.step(delta).0)
    }
    /// `Bᵀ Log(E)`, `E` the minimal rotation taking `n` to `other` (CH.9(b)).
    ///
    /// # Domain
    ///
    /// `other ≠ -n`, where every `δ` of norm `π` is a preimage and `local` has no limit (CH.9(d));
    /// on the far hemisphere (`n·m <= 0`), `‖n × m‖²` normal, so `m` farther from `-n` than about
    /// `1.5e-154` (`f64`) or `1.1e-19` (`f32`). Nothing is asserted: at `m = -n` the answer is
    /// `0`, and closer than the bound a finite value of the wrong size (`NUMERICS.md` §12).
    #[inline]
    fn local(&self, other: &S2<S>) -> RnTangent<S, 2> {
        self.log_parts(other).delta
    }
    /// `B'ᵀ J_l(Bδ) B` (CH.10(a)), `B'` the basis at `m = retract(δ)`, the same bits `retract`
    /// returns, so the frame it reads `m` in is `S2Chart::at(&retract(δ))`'s.
    ///
    /// By CH.8(c) (`b₁ × b₂ = ς n`), `J_l b₁ = (1 - bδ₂²) b₁ + bδ₁δ₂ b₂ - aςδ₂ n` and
    /// `J_l b₂ = bδ₁δ₂ b₁ + (1 - bδ₁²) b₂ + aςδ₁ n`: the symmetric `A(δ)` of CH.10 plus the normal
    /// part, which `B'` reads and `B` does not. Defined for every `δ`; singular at `‖δ‖ = π`, where
    /// every direction retracts to `-n` (CH.10(a)).
    #[inline]
    fn retract_jacobian(&self, delta: &RnTangent<S, 2>) -> Mat2<S> {
        let (m, z) = self.step(delta);
        let (a, b) = jr_coeffs(z);
        let [d1, d2] = delta.rho.0;
        let n = self.base.0;
        let e = sign(n) * a;
        let bx = b * d1 * d2;
        let col = |x: S, y: S, t: S| {
            Vector(array::from_fn(|k| {
                (x * self.b1.0[k] + y * self.b2.0[k]) + t * n.0[k]
            }))
        };
        let j1 = col(S::one() - b * d2 * d2, bx, -(e * d2));
        let j2 = col(bx, S::one() - b * d1 * d1, e * d1);
        let (c1, c2) = basis(m);
        Matrix::from_rows([
            Vector([c1.dot(j1), c1.dot(j2)]),
            Vector([c2.dot(j1), c2.dot(j2)]),
        ])
    }
    /// `A(δ)⁻¹ Z` at `δ = local(other)` (CH.10(b)), the inverse of `retract_jacobian(δ)`.
    ///
    /// # Domain
    ///
    /// [`local`](Chart::local)'s. Its larger singular value is `α/sin α`, which grows without
    /// bound as `other → -n` (CH.10(b)).
    #[inline]
    fn local_jacobian(&self, other: &S2<S>) -> Mat2<S> {
        let l = self.log_parts(other);
        self.local_jac(other, &l)
    }
    /// One `log_ratio` and one cross product for both (`0060` decision 3).
    ///
    /// # Domain
    ///
    /// [`local`](Chart::local)'s.
    #[inline]
    fn local_with_jacobian(&self, other: &S2<S>) -> (RnTangent<S, 2>, Mat2<S>) {
        let l = self.log_parts(other);
        (l.delta, self.local_jac(other, &l))
    }
}

/// S² linearizes in its Householder chart; its coordinates are an [`RnTangent<S, 2>`]
/// (`0067`).
impl<S: Real> Manifold<S> for S2<S> {
    type Tangent = RnTangent<S, 2>;
    type Chart = S2Chart<S>;
    const DOF: usize = 2;
}

/// `B rows · v`, the two rows of `m` against the first two entries of `v`.
#[inline]
fn apply_rows<S: Real, const D: usize>(m: &Mat2<S>, v: &Vector<S, D>) -> Vector<S, D> {
    let rows = [m.row(0).0, m.row(1).0];
    Vector(array::from_fn(|r| {
        rows[r][0] * v.0[0] + rows[r][1] * v.0[1]
    }))
}

/// The dense Jacobian of [`S2Chart`] (`docs/PHASE5.md` §1.4): S² has no group structure for a
/// structured type to exploit.
impl<S: Real> Jac<S, RnTangent<S, 2>> for Mat2<S> {
    #[inline]
    fn identity() -> Self {
        Matrix::identity()
    }
    #[inline]
    fn mul(&self, o: &Self) -> Self {
        *self * *o
    }
    /// The algebraic inverse, `Mat2::inverse_adj`: the adjugate, exact, over the determinant
    /// (`0005`, `0061`).
    ///
    /// # Domain
    ///
    /// `det != 0`, checked by `debug_assert!` ("no lane has `det = 0`", so a NaN matrix returns
    /// NaN); `retract_jacobian` at `‖δ‖ = π` is singular. A release build divides and returns
    /// `±inf` or NaN (D11).
    #[inline]
    fn inverse(&self) -> Self {
        let (inv, det) = self.inverse_adj();
        debug_assert!(
            !det.abs().le(S::zero()).any(),
            "Jac::inverse: singular 2 x 2"
        );
        inv
    }
    #[inline]
    fn neg(&self) -> Self {
        -*self
    }
    #[inline]
    fn apply(&self, t: &RnTangent<S, 2>) -> RnTangent<S, 2> {
        RnTangent { rho: *self * t.rho }
    }
    #[inline]
    fn apply_transpose(&self, t: &RnTangent<S, 2>) -> RnTangent<S, 2> {
        RnTangent {
            rho: self.transpose() * t.rho,
        }
    }
    #[inline]
    fn write_dense(&self, out: &mut StridedMut<'_, S>) {
        debug_assert!(
            out.rows() == 2 && out.cols() == 2,
            "Jac::write_dense: the view is not DOF x DOF"
        );
        for c in 0..2 {
            for (r, v) in self.col(c).0.iter().enumerate() {
                out.set(r, c, *v);
            }
        }
    }
    #[inline]
    fn sandwich<const D: usize>(&self, cov: &Matrix<S, D, D>) -> Matrix<S, D, D> {
        const { assert!(D == <RnTangent<S, 2> as Tangent<S>>::DOF) };
        // `Mat3`'s order: `M = J Σ` by columns, then row `r` of `M Jᵀ` is `J` applied to row `r`.
        let m = Matrix::from_cols(array::from_fn(|c| apply_rows::<S, D>(self, &cov.col(c))));
        Matrix::from_rows(array::from_fn(|r| apply_rows::<S, D>(self, &m.row(r))))
    }
}
