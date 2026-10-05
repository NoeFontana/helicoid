//! SE_N(3): the group `SEn3`, its tangent `SEn3Tangent` and the twist converters
//! (`docs/PHASE3.md` §5). The structured Jacobian `SEn3Jac` is `dualmat`'s.

use crate::coeffs::{jr_coeffs, jr_inv_coeff, q_coeffs};
use crate::dualmat::SEn3Jac;
use crate::quat::Quat;
use crate::side::Side;
use crate::so3::{norm_sq, SO3Tangent, SO3};
use crate::traits::{tie_dof, Jac, LieGroup, Tangent};
use core::array;
use core::iter::once;
use core::ops::Mul;
use helicoid_linalg::{hat, Blend, Mat3, Matrix, Point, Point3, Real, Vec3, Vector};

/// A tangent vector of SE_N(3): `[φ; ρ₁; …; ρ_N]`, rotation first (`NUMERICS.md` §1, `0002`).
///
/// The fields are named (`docs/API.md` R3). The flat order, `phi` and then `rho[0]`, `rho[1]`, …,
/// is visible only through [`write_dense`](Tangent::write_dense) and
/// [`read_dense`](Tangent::read_dense) and, for a [`Twist`], the two converters. `DOF` is
/// `3 + 3N`. The operations are the methods of [`Tangent`]: there is no `Add`, `Sub`, `Index` or
/// `From<[S; 6]>` (R1, R3), and none of these compiles:
///
/// ```compile_fail,E0369
/// use helicoid::Twist;
/// use helicoid_linalg::Vector;
/// let t = Twist { phi: Vector([1.0_f64; 3]), rho: [Vector([1.0; 3])] };
/// let _ = t + t;
/// ```
///
/// ```compile_fail,E0608
/// use helicoid::Twist;
/// use helicoid_linalg::Vector;
/// let t = Twist { phi: Vector([1.0_f64; 3]), rho: [Vector([1.0; 3])] };
/// let _ = t[0];
/// ```
///
/// ```compile_fail,E0277
/// use helicoid::Twist;
/// let _: Twist<f64> = [0.0_f64; 6].into();
/// ```
///
/// Positive control, so the failures above come from the missing impls and not from an import:
///
/// ```
/// use helicoid::{Tangent, Twist};
/// use helicoid_linalg::Vector;
/// let t = Twist { phi: Vector([1.0_f64; 3]), rho: [Vector([1.0; 3])] };
/// let mut out = [0.0; 6];
/// t.add(&t).write_dense(&mut out);
/// assert_eq!(out, [2.0; 6]);
/// ```
#[derive(Clone, Copy, Debug)]
pub struct SEn3Tangent<S, const N: usize> {
    /// The rotation part `φ`; `θ = ‖φ‖`.
    pub phi: Vec3<S>,
    /// The translation parts `ρ₁, …, ρ_N`, one per column of the group element. For SE₂(3) they
    /// are the velocity part and then the position part (`docs/maths/se3.md` SE.1).
    pub rho: [Vec3<S>; N],
}

/// The tangent of SE(3): `omega()` is `phi` and `v()` is `rho[0]`.
pub type Twist<S> = SEn3Tangent<S, 1>;

impl<S: Real, const N: usize> SEn3Tangent<S, N> {
    /// The `3 + 3N` components in the dense order.
    fn comps(self) -> impl Iterator<Item = S> {
        once(self.phi).chain(self.rho).flat_map(|v| v.0)
    }

    /// The tangent whose dense component `i` is `at(i)`, in the order [`comps`](Self::comps) gives.
    ///
    /// One body for both arms of [`read_dense`](Tangent::read_dense), so the in-range arm and the
    /// NaN-poisoning one cannot drift apart in the order they assign.
    #[inline]
    fn from_dense_fn(at: impl Fn(usize) -> S) -> Self {
        Self {
            phi: Vector(array::from_fn(&at)),
            rho: array::from_fn(|i| Vector(array::from_fn(|r| at(3 * (i + 1) + r)))),
        }
    }
}

impl<S: Real, const N: usize> Blend<S> for SEn3Tangent<S, N> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self {
            phi: Vector::blend(m, t.phi, f.phi),
            rho: <[Vec3<S>; N]>::blend(m, t.rho, f.rho),
        }
    }
}

impl<S: Real, const N: usize> Tangent<S> for SEn3Tangent<S, N> {
    const DOF: usize = 3 + 3 * N;
    #[inline]
    fn zero() -> Self {
        Self {
            phi: Vector([S::zero(); 3]),
            rho: [Vector([S::zero(); 3]); N],
        }
    }
    #[inline]
    fn add(&self, o: &Self) -> Self {
        Self {
            phi: self.phi + o.phi,
            rho: array::from_fn(|i| self.rho[i] + o.rho[i]),
        }
    }
    #[inline]
    fn sub(&self, o: &Self) -> Self {
        Self {
            phi: self.phi - o.phi,
            rho: array::from_fn(|i| self.rho[i] - o.rho[i]),
        }
    }
    #[inline]
    fn neg(&self) -> Self {
        Self {
            phi: -self.phi,
            rho: self.rho.map(|v| -v),
        }
    }
    #[inline]
    fn scale(&self, k: S) -> Self {
        Self {
            phi: self.phi.scale(k),
            rho: self.rho.map(|v| v.scale(k)),
        }
    }
    // One left-to-right sum over all `3 + 3N` products, not a sum of per-block dots (D16). The
    // accumulator, not `dot`, is the required operation (`0025`), and the `+0` the provided `dot`
    // seeds it with is normative: seeding from the first product instead differs on a signed zero.
    #[inline]
    fn dot_acc(&self, o: &Self, acc: S) -> S {
        // Indexed, not over `comps`: the order is the one `comps` gives and the products are the
        // same, with no `Chain`/`FlatMap` state for the accumulator to carry.
        let mut a = acc;
        for r in 0..3 {
            a = a + self.phi.0[r] * o.phi.0[r];
        }
        for i in 0..N {
            for r in 0..3 {
                a = a + self.rho[i].0[r] * o.rho[i].0[r];
            }
        }
        a
    }
    #[inline]
    fn write_dense(&self, out: &mut [S]) {
        debug_assert!(out.len() == Self::DOF, "Tangent::write_dense: wrong length");
        // Indexed, as `dot_acc` is and for its reason: the order is the one `comps` gives, with no
        // `Chain`/`FlatMap` state to carry on the path every dense export and every `apply_flat`
        // column takes. The one `get_mut` bounds the whole arm, which then stores at constant
        // a slice of length `DOF`, which `copy_from_slice` then fills per block: at `N = 3` that arm
        // is one `cmp`/`jb` and a `movups`-packed copy in release on x86_64, with no bounds branch
        // and no `memcpy` call, against twelve scalar stores driven by a `Chain<Once, IntoIter>`
        // inside a `FlatMap`.
        // A view longer or shorter than `DOF` keeps what the `zip` did -- the prefix that fits, the
        // tail untouched -- which `out_of_domain_does_not_panic_in_release` pins.
        if let Some(d) = out.get_mut(..Self::DOF) {
            d[..3].copy_from_slice(&self.phi.0);
            for i in 0..N {
                d[3 * (i + 1)..3 * (i + 2)].copy_from_slice(&self.rho[i].0);
            }
        } else {
            for (o, v) in out.iter_mut().zip(self.comps()) {
                *o = v;
            }
        }
    }
    #[inline]
    fn read_dense(src: &[S]) -> Self {
        debug_assert!(src.len() == Self::DOF, "Tangent::read_dense: wrong length");
        // A short `src` is out of domain and must not produce a usable tangent: `+0` is a valid
        // component, so it would hand a solver a plausible wrong update, while NaN propagates to
        // whatever the caller computes. D11 forbids the release check that would say so instead.
        //
        // One length test, not one per component. `ProductJac::sandwich` reaches `read_dense` once
        // per column and once per row of its `D x D` argument, and again per nesting level, so the
        // in-range arm reads a slice of length `DOF` at constant indices and carries no per-entry
        // bound: at `N = 3` it is one `cmp`/`jb` and twelve loads in release on x86_64. The
        // poisoning arm below is the shape a per-component `get` had on every call -- thirteen
        // compares and as many selects -- and it is now only reached out of domain. Only the
        // missing entries are poisoned, which is `0025` decision 4 and what
        // `out_of_domain_does_not_panic_in_release` pins.
        match src.get(..Self::DOF) {
            Some(d) => Self::from_dense_fn(|i| d[i]),
            None => Self::from_dense_fn(|i| {
                src.get(i).copied().unwrap_or_else(|| S::zero() / S::zero())
            }),
        }
    }
}

impl<S: Real> SEn3Tangent<S, 1> {
    /// The angular velocity `ω`, the rotation part.
    #[inline]
    pub fn omega(&self) -> Vec3<S> {
        self.phi
    }

    /// The linear velocity `v`, the translation part.
    #[inline]
    pub fn v(&self) -> Vec3<S> {
        self.rho[0]
    }

    /// The twist of a translation-first array `[v; ω]`, the order of Barfoot (`[ρ; φ]`), Solà,
    /// Sophus and manif (`[υ; ω]`), and of locus-tag's `delta` (`0002`).
    ///
    /// It is the permutation `Π = [[0, I₃], [I₃, 0]]` of `docs/maths/se3.md` SE.14(a), `[v; ω]`
    /// to `[ω; v]`: `a[0..3]` becomes `v` and `a[3..6]` becomes `ω`. Six scalars are copied, so
    /// the conversion is exact. It converts the order of a tangent and nothing else: a covariance
    /// or a Jacobian in the other order is `Π M Πᵀ`, the caller's (SE.14(b)), and locus-tag's
    /// decoupled retraction is a different map from `Exp`, not `Π` of it (`0012`).
    ///
    /// ```
    /// use helicoid::Twist;
    /// let t = Twist::from_translation_first([1.0_f64, 2.0, 3.0, 4.0, 5.0, 6.0]);
    /// assert_eq!((t.v().0, t.omega().0), ([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]));
    /// ```
    #[inline]
    pub fn from_translation_first(a: [S; 6]) -> Self {
        Self {
            phi: Vector([a[3], a[4], a[5]]),
            rho: [Vector([a[0], a[1], a[2]])],
        }
    }

    /// The translation-first array `[v; ω]`: the inverse of
    /// [`from_translation_first`](SEn3Tangent::from_translation_first), and the same
    /// permutation `Π` (which is its own inverse for `N = 1`). Exact.
    #[inline]
    pub fn to_translation_first(&self) -> [S; 6] {
        let ([w0, w1, w2], [v0, v1, v2]) = (self.phi.0, self.rho[0].0);
        [v0, v1, v2, w0, w1, w2]
    }
}

/// An element of SE_N(3): a rotation and `N` translation columns, `(R, x₁, …, x_N)`
/// (`docs/PHASE3.md` §5).
///
/// `N = 1` is SE(3), a rigid transform ([`SE3`]); `N = 2` is SE₂(3), the InEKF extended pose
/// ([`SE23`]), whose columns are the velocity and then the position (`docs/maths/se3.md` SE.1).
/// Composition is `a * b = T_a_x · T_x_b` (`0002`) and the action is `Mul<Point3<S>>` at `N = 1`.
///
/// The fields are private and `repr(C)` is `PHASE3.md` §5's spelling, not a promise: layout is not
/// a semver contract (D2) and a consumer keeps its own storage format. The stored quaternion
/// carries the unit invariant, guarded where a raw one enters — [`SE3::from_quat_translation`]'s
/// `debug_assert!` — and nowhere else, since `Exp` and composition produce it. A value is built by
/// [`exp`](LieGroup::exp), [`SE3::from_rt`] or [`SE3::from_quat_translation`] and read by
/// [`rotation`](Self::rotation), the `N`-specific accessors and [`log`](LieGroup::log).
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct SEn3<S, const N: usize> {
    q: Quat<S>,
    x: [Vec3<S>; N],
}

/// SE(3), a rigid transform: `N = 1`, whose tangent is a [`Twist`].
pub type SE3<S> = SEn3<S, 1>;

/// SE₂(3), the InEKF extended pose: `N = 2`, the columns `x₁ = v` then `x₂ = p`.
pub type SE23<S> = SEn3<S, 2>;

impl<S: Real, const N: usize> SEn3<S, N> {
    /// The rotation part `R`.
    #[inline]
    pub fn rotation(&self) -> SO3<S> {
        SO3::from_quat_unchecked(self.q)
    }

    /// The element of a rotation and its columns, exactly: a move of `3 + 4` numbers, no
    /// transcendental and no round trip through `Exp`/`Log` (`0042` draft, which is `0029` one
    /// level down — the same question `Product::from_parts` answered, for the concrete group).
    ///
    /// The named forms are the ones to prefer where they exist: [`SE3::from_rt`],
    /// [`SE3::from_quat_translation`]. This one is what a general `N` has, and what the
    /// conformance subject reads the `sen3_*_n3` corpus records with.
    #[inline]
    pub fn from_parts(r: SO3<S>, x: [Vec3<S>; N]) -> Self {
        Self { q: r.quat(), x }
    }

    /// The rotation and the columns, exactly: the inverse of [`from_parts`](Self::from_parts), and
    /// the only way to read an `N > 2` element without a `Log` (`0042` draft).
    #[inline]
    pub fn parts(&self) -> (SO3<S>, [Vec3<S>; N]) {
        (self.rotation(), self.x)
    }
}

impl<S: Real> SEn3<S, 1> {
    /// The translation `t`.
    #[inline]
    pub fn translation(&self) -> Vec3<S> {
        self.x[0]
    }

    /// The rigid transform `(R, t)`.
    #[inline]
    pub fn from_rt(r: SO3<S>, t: Vec3<S>) -> Self {
        Self::from_parts(r, [t])
    }

    /// The rigid transform of a quaternion and a translation.
    ///
    /// The normalizing path is `from_rt` of [`SO3::from_quat_normalized`], which divides by the
    /// norm (`0027`); this one does not.
    ///
    /// # Domain
    ///
    /// `q` is unit within `NUMERICS.md` §3.6's bound, by `debug_assert!` and no release check
    /// (D11). The assert is reached through [`Quat::from_wxyz_unchecked`], which is the item that
    /// carries it — [`SO3::from_quat_unchecked`] is a move and checks nothing, so routing through
    /// it would have made this constructor's stated domain unenforced.
    #[inline]
    pub fn from_quat_translation(q: Quat<S>, t: Vec3<S>) -> Self {
        let q = Quat::from_wxyz_unchecked(q.w, q.x, q.y, q.z);
        Self::from_rt(SO3::from_quat_unchecked(q), t)
    }

    /// `X p = R p + t` (`NUMERICS.md` §5.5), the rotation through the quaternion sandwich §3.3
    /// states, as [`SO3::act`] forms it.
    #[inline]
    fn act(&self, p: Point3<S>) -> Point3<S> {
        Point((self.rotation().act(Vector(p.0)) + self.x[0]).0)
    }

    /// `X p` for every point in `pts`, in place: one `R(q)` for the whole slice, as
    /// [`SO3::act_many`] forms it, then one translation each.
    ///
    /// The matrix rounds its nine entries before any point is touched, where the `Mul<Point3>`
    /// action rounds a sandwich per point, so the two differ in the bits — the difference §14's
    /// `act_many` row records for SO(3).
    #[inline]
    pub fn act_many(&self, pts: &mut [Point3<S>]) {
        let r = self.rotation().to_matrix();
        let t = self.x[0];
        for p in pts.iter_mut() {
            *p = Point((r * Vector(p.0) + t).0);
        }
    }

    /// `(∂(X p)/∂X, ∂(X p)/∂p)` in the convention of `Sd` (`NUMERICS.md` §2.4): rotation-first
    /// `[−R [p]_×, R]` and `R` on the right, `[−[X p]_×, I]` and `R` on the left.
    ///
    /// The `3 x 6` is the one place a group method returns a dense matrix (`0005` permits it: it
    /// is not a `Jac`, the tangent of no group is its row space, and §2.4 states it dense).
    #[inline]
    pub fn act_jacobians<Sd: Side>(&self, p: Point3<S>) -> (Matrix<S, 3, 6>, Mat3<S>) {
        let r = self.rotation().to_matrix();
        let (left, right) = match Sd::IS_RIGHT {
            true => (-(r * hat(Vector(p.0))), r),
            false => (-hat(Vector(self.act(p).0)), Mat3::identity()),
        };
        let cols = array::from_fn(|c| match c < 3 {
            true => left.col(c),
            false => right.col(c - 3),
        });
        (Matrix::from_cols(cols), r)
    }
}

impl<S: Real> SEn3<S, 2> {
    /// The velocity `v`, the first column (`docs/maths/se3.md` SE.1).
    #[inline]
    pub fn velocity(&self) -> Vec3<S> {
        self.x[0]
    }

    /// The position `p`, the second column (`docs/maths/se3.md` SE.1).
    #[inline]
    pub fn position(&self) -> Vec3<S> {
        self.x[1]
    }
}

impl<S: Real, const N: usize> Blend<S> for SEn3<S, N> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self {
            q: Quat::blend(m, t.q, f.q),
            x: <[Vec3<S>; N]>::blend(m, t.x, f.x),
        }
    }
}

impl<S: Real, const N: usize> Mul for SEn3<S, N> {
    type Output = Self;

    /// `a · b = (R_a R_b, x_{a,i} + R_a x_{b,i})`: the Hamilton product of the two quaternions,
    /// and the sandwich of §3.3 applied to each of `b`'s columns.
    ///
    /// The product quaternion is not renormalized — `Exp` and composition both return what the
    /// arithmetic gives, and `SO3::renormalize` is the caller's step (`0027`).
    #[inline]
    // The `+` inside a `Mul` is the group law, not a slip: a rigid motion composes by rotating the
    // second element's columns and *adding* the first's (`NUMERICS.md` §5, `0002`).
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn mul(self, o: Self) -> Self {
        let r = self.rotation();
        Self {
            q: self.q * o.q,
            x: array::from_fn(|i| self.x[i] + r.act(o.x[i])),
        }
    }
}

impl<S: Real> Mul<Point3<S>> for SEn3<S, 1> {
    type Output = Point3<S>;

    /// `X p = R p + t` (`NUMERICS.md` §5.5). The operator is the action (`docs/PHASE3.md` §5): a
    /// [`Point3`] is acted on, a `Vec3` is not, which is what keeps the translation out of the
    /// rotation of a direction.
    #[inline]
    fn mul(self, p: Point3<S>) -> Point3<S> {
        self.act(p)
    }
}

/// Barfoot's block `Q(ρ, φ)` of `NUMERICS.md` §5.3, with `W = φ^` and `X = ρ^` given as matrices
/// and `(b, d, e)` from one grouped `q_coeffs` branch.
///
/// Every word is associated as §5.3 writes it and the products are the seeded kernel's, operand
/// for operand (`docs/maths/coefficients.md` CO.6, `0010`): the `xtask` subject that scores the
/// `sen3_j*` corpus ids runs this same sequence, so a parity gap between them is a coefficient's
/// and not an association's.
///
/// `W²` is passed in, not formed here: it is the one word of the eight that does not depend on
/// `ρ`, so forming it per column cost `(N − 1) · 27` multiplications for the same bits.
#[inline]
fn q_block<S: Real>(x: &Mat3<S>, w: &Mat3<S>, ww: &Mat3<S>, b: S, d: S, e: S) -> Mat3<S> {
    let (wx, xw, ww) = (*w * *x, *x * *w, *ww);
    let (wxw, wwx, xww) = (wx * *w, ww * *x, xw * *w);
    let b_words = (wx + xw) + wxw;
    let d_words = (wwx + xww) - wxw.scale(S::lit(3.0));
    let e_words = wxw * *w + wwx * *w;
    let head = x.scale(S::lit(0.5)) + b_words.scale(b);
    (head + d_words.scale(d)) + e_words.scale(e)
}

impl<S: Real, const N: usize> LieGroup<S> for SEn3<S, N> {
    type Tangent = SEn3Tangent<S, N>;
    type Jac = SEn3Jac<S, N>;
    const DOF: usize = 3 + 3 * N;

    #[inline]
    fn identity() -> Self {
        const { tie_dof::<S, Self>() };
        Self::from_parts(SO3::identity(), [Vector([S::zero(); 3]); N])
    }

    /// `X⁻¹ = (R⁻¹, −R⁻¹ x_i)`: the conjugate quaternion, and its sandwich applied to each
    /// column, negated.
    #[inline]
    fn inverse(&self) -> Self {
        let r = self.rotation().inverse();
        Self::from_parts(r, self.x.map(|v| -r.act(v)))
    }

    /// `Exp(τ) = (Exp(φ), J_l(φ) ρ_i)` (`NUMERICS.md` §5.1).
    ///
    /// `J_l(φ) ρ = ρ + a (φ × ρ) + b (φ × (φ × ρ))`, two cross products per column rather than the
    /// assembled `I + aW + bW²` times `ρ`. The two are the same map and different roundings, and
    /// `0036` (draft) measured the choice on the real corpus: 16 domination failures against 36 for
    /// the matrix form on `sen3_exp_n1`, the cross form better away from `π`. It is also the
    /// cheaper one at every `N` this type has (`12N` multiplications against `27 + 9N`) and it is
    /// what `tf_tree_math::exp_se3` does, so the migration's parity rows compare one program's
    /// rounding with its own (`0010`, `0041` draft).
    #[inline]
    fn exp(tau: &SEn3Tangent<S, N>) -> Self {
        let phi = tau.phi;
        let (a, b) = jr_coeffs(norm_sq(phi));
        Self {
            q: SO3::exp(&SO3Tangent { phi }).quat(),
            x: tau.rho.map(|p| {
                let c1 = phi.cross(p);
                (p + c1.scale(a)) + phi.cross(c1).scale(b)
            }),
        }
    }

    /// `Log(X) = (φ = Log(R), ρ_i = J_l⁻¹(φ) x_i)` (`NUMERICS.md` §5.1).
    ///
    /// `J_l⁻¹(φ) x = x − ½ (φ × x) + c (φ × (φ × x))`, the `J_l⁻¹ = J_r⁻¹(−φ)` of §3.5 applied by
    /// two cross products, for [`exp`](LieGroup::exp)'s reasons and matching
    /// `tf_tree_math::log_se3`.
    ///
    /// # Domain
    ///
    /// `θ < 2π`, as [`SO3::jr_inv`]'s is: `c` has a pole there. Nothing is asserted, for that
    /// method's reason.
    #[inline]
    fn log(&self) -> SEn3Tangent<S, N> {
        let phi = self.rotation().log().phi;
        let c = jr_inv_coeff(norm_sq(phi));
        let half = S::lit(0.5);
        SEn3Tangent {
            phi,
            rho: self.x.map(|v| {
                let c1 = phi.cross(v);
                (v - c1.scale(half)) + phi.cross(c1).scale(c)
            }),
        }
    }

    /// `Ad_X = R + ε [x_i]_× R` (`NUMERICS.md` §5.2).
    #[inline]
    fn adjoint(&self) -> SEn3Jac<S, N> {
        let r = self.rotation().to_matrix();
        SEn3Jac {
            diag: r,
            col: self.x.map(|v| hat(v) * r),
        }
    }

    /// `ad_τ = W + ε [ρ_i]_×` (`NUMERICS.md` §5.2).
    #[inline]
    fn ad(tau: &SEn3Tangent<S, N>) -> SEn3Jac<S, N> {
        SEn3Jac {
            diag: hat(tau.phi),
            col: tau.rho.map(hat),
        }
    }

    /// `J_r(τ) = J_r(φ) + ε Q(−ρ_i, −φ)` (`NUMERICS.md` §5.3), the diagonal block from
    /// [`SO3::jr`] and the columns from this module's `q_block`.
    ///
    /// `jl` is **not** overridden: the trait's provided `J_r(−τ)` is `J_l(φ) + ε Q(ρ_i, φ)`,
    /// which is §5.3's left form exactly — negating the tangent negates both arguments of `Q`
    /// twice over — so one code path serves both sides, as it does for [`SO3`].
    #[inline]
    fn jr(tau: &SEn3Tangent<S, N>) -> SEn3Jac<S, N> {
        let (b, d, e) = q_coeffs(norm_sq(tau.phi));
        let w = hat(-tau.phi);
        let ww = w * w;
        SEn3Jac {
            diag: SO3::jr(&SO3Tangent { phi: tau.phi }),
            col: tau.rho.map(|r| q_block(&hat(-r), &w, &ww, b, d, e)),
        }
    }

    /// `J_r⁻¹(τ) = J_r⁻¹(φ) − ε J_r⁻¹(φ) Q(−ρ_i, −φ) J_r⁻¹(φ)` (`NUMERICS.md` §5.4).
    ///
    /// The dual-matrix inverse of [`jr`](LieGroup::jr), with `A⁻¹` taken from §3.5's **closed
    /// form** — [`SO3::jr_inv`] — and not from [`Jac::inverse`], which would reach
    /// `Mat3::inverse_adj` and divide by a computed `det A`. That is what §5.4 states, and it is
    /// what lets this method keep §12's promise of a finite value: the closed form has no
    /// division at all, so there is no `det A = 0` for a release build to divide by. `Jac::inverse`
    /// remains the inverse of a `SEn3Jac` that is nobody's `J_r` (`0005`).
    ///
    /// # Domain
    ///
    /// `θ < 2π`, [`SO3::jr_inv`]'s, asserted nowhere for its reason.
    #[inline]
    fn jr_inv(tau: &SEn3Tangent<S, N>) -> SEn3Jac<S, N> {
        let ai = SO3::jr_inv(&SO3Tangent { phi: tau.phi });
        let (b, d, e) = q_coeffs(norm_sq(tau.phi));
        let w = hat(-tau.phi);
        let ww = w * w;
        SEn3Jac {
            diag: ai,
            col: tau
                .rho
                .map(|r| -(ai * q_block(&hat(-r), &w, &ww, b, d, e) * ai)),
        }
    }

    /// `(Ad_Exp(τ)⁻¹, J_r(τ))` (`NUMERICS.md` §2.3).
    ///
    /// `Ad_Exp(τ)⁻¹` as `Ad_Exp(−τ)`, not as [`Jac::inverse`] of `Ad` and not as `Ad` of the group
    /// inverse. All three are the same dual matrix — `Ad` is a homomorphism and `Exp(−τ)` is
    /// `Exp(τ)⁻¹` — and this one is the cheapest *and* the exact one: `Exp(−τ)`'s quaternion is
    /// `Exp(τ)`'s conjugate to the bit (the coefficients are even in `θ` and a negation is exact),
    /// so the saving is the group inverse's `N` sandwiches, and the result is bit-identical to the
    /// row `laws::jacobian_rows` checks against, which is why `rows` is recorded at exactly `0`.
    /// Through `Ad` of the inverse it was 2.182 `u`: `inverse` reaches its columns as `−Rᵗ(J_l ρ)`
    /// where this reaches them as `−J_r ρ`, equal in exact arithmetic and not in the bits.
    #[inline]
    fn rplus_jacobians(&self, tau: &SEn3Tangent<S, N>) -> (SEn3Jac<S, N>, SEn3Jac<S, N>) {
        (Self::exp(&tau.neg()).adjoint(), Self::jr(tau))
    }

    /// `(Ad_Exp(τ), J_l(τ))` (`NUMERICS.md` §2.3).
    #[inline]
    fn lplus_jacobians(&self, tau: &SEn3Tangent<S, N>) -> (SEn3Jac<S, N>, SEn3Jac<S, N>) {
        (Self::exp(tau).adjoint(), Self::jl(tau))
    }

    /// `(J_r⁻¹(τ), −J_l⁻¹(τ))` at `τ = self ⊖_R base` (`NUMERICS.md` §2.3).
    ///
    /// Two inversions, where [`SO3::rminus_jacobians`] takes one and transposes it. The *diagonal*
    /// blocks do transpose into each other exactly, as that method's note proves, but the columns
    /// do not: `J_l⁻¹(τ)`'s is `Q(ρ, φ)`'s and `J_r⁻¹(τ)`'s is `Q(−ρ, −φ)`'s, and the transpose of
    /// a `SEn3Jac` is block *upper* triangular, so it is not a `SEn3Jac` to return.
    ///
    /// A fused body is available and not taken here: the two share `θ²`, both coefficient
    /// branches, `W²`, and six of `Q`'s eight matrix words up to an exact sign, so roughly half of
    /// this call is recomputation on what `SO3::rminus_jacobians` calls the crate's hottest
    /// Jacobian path. It is a measurement `PHASE3.md` §11's benches are owed — and a second
    /// rounding to record against §14 — not a reading of §2.3, so it waits for the bench.
    #[inline]
    fn rminus_jacobians(&self, base: &Self) -> (SEn3Jac<S, N>, SEn3Jac<S, N>) {
        let tau = self.rminus(base);
        (Self::jr_inv(&tau), Self::jl_inv(&tau).neg())
    }

    /// `(J_l⁻¹(τ), −J_r⁻¹(τ))` at `τ = self ⊖_L base` (`NUMERICS.md` §2.3); see
    /// [`rminus_jacobians`](LieGroup::rminus_jacobians) for the two inversions.
    #[inline]
    fn lminus_jacobians(&self, base: &Self) -> (SEn3Jac<S, N>, SEn3Jac<S, N>) {
        let tau = self.lminus(base);
        (Self::jl_inv(&tau), Self::jr_inv(&tau).neg())
    }

    /// Right `(Ad_Y⁻¹, I)`, left `(I, Ad_X)` (`NUMERICS.md` §2.3), selected by `Sd::IS_RIGHT` at
    /// monomorphization; `Ad⁻¹` as in [`rplus_jacobians`](LieGroup::rplus_jacobians).
    #[inline]
    fn compose_jacobians<Sd: Side>(&self, rhs: &Self) -> (SEn3Jac<S, N>, SEn3Jac<S, N>) {
        match Sd::IS_RIGHT {
            true => (rhs.inverse().adjoint(), SEn3Jac::identity()),
            false => (SEn3Jac::identity(), self.adjoint()),
        }
    }

    /// Right `−Ad_X`, left `−Ad_X⁻¹` (`NUMERICS.md` §2.3).
    #[inline]
    fn inverse_jacobian<Sd: Side>(&self) -> SEn3Jac<S, N> {
        match Sd::IS_RIGHT {
            true => self.adjoint().neg(),
            false => self.inverse().adjoint().neg(),
        }
    }
}
