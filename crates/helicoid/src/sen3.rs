//! SE_N(3): the group `SEn3`, its tangent `SEn3Tangent` and the twist converters
//! (`docs/PHASE3.md` §5). The structured Jacobian `SEn3Jac` is `dualmat`'s.

use crate::coeffs::log_ratio_short;
use crate::coeffs::{jr_coeffs, jr_coeffs_sinc, jr_inv_coeff, jr_inv_coeff_dz, q_five_coeffs};
use crate::dualmat::{zero3, SEn3Jac};
use crate::quat::Quat;
use crate::side::Side;
use crate::so3::{
    geodesic_long, geodesic_rel, geodesic_short, hat_mul, mul_hat, norm_sq, HatSq, SO3Tangent, SO3,
};
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
    /// Sophus, and manif (`[υ; ω]`) ([`0002`](./docs/decisions/0002-one-convention-for-a-stack-that-already-disagrees.md)).
    ///
    /// Applies permutation `Π = [[0, I₃], [I₃, 0]]` (`docs/maths/se3.md` SE.14(a)) from `[v; ω]`
    /// to `[ω; v]`. Covariances or Jacobians in that order transform as `Π M Πᵀ` (SE.14(b)).
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
    /// [`Quat::renormalize`] on the rotation's quaternion, in place, leaving the columns alone:
    /// the step `Mul`'s rustdoc names as the caller's, in one call rather than a
    /// `parts`/`from_parts` round trip around a private field (`0044` item 2).
    ///
    /// # Domain
    ///
    /// None, as [`Quat::renormalize`]'s is none: the step is defined for every input and only its
    /// *accuracy* has a domain — one step is a normalization for `|‖q‖² − 1|` at most `2^-26.29`
    /// (`f64`) or `2^-11.79` (`f32`), the drift a chain of compositions accumulates.
    #[inline]
    pub fn renormalize(&mut self) {
        self.q.renormalize();
    }

    /// `self · other⁻¹`, without forming `other⁻¹` (`0048`).
    ///
    /// `(q_a q_b*, x_{a,i} − R(q_a q_b*) x_{b,i})`: **one** rotation per column, where
    /// `*self * other.inverse()` does two — `inverse` forms `−R_bᵗ x_{b,i}` and the product then
    /// rotates it by `R_a` — because `R_a R_bᵗ` is `R(q_a q_b*)`. It is [`lminus`]'s first half,
    /// which is why that method's body is this one (`0048` decision 3).
    ///
    /// The two are equal in exact arithmetic and **not bit-identical**, so this is its own routine
    /// with its own `NUMERICS.md` §14 twin, whose reference *is* the composition it differs from
    /// (`sen3_tests::group::sen3_mul_inv_matches_reference_n1`).
    ///
    /// The agreement is `O(u)` **against the input scale**, not against the result: both forms
    /// carry an absolute error of order `u·max‖x‖`, so where the two frames nearly coincide and
    /// the difference cancels, the *relative* gap grows as that cancellation. Measured worst
    /// `8.9 u` of `max(‖x_a‖, ‖x_b‖, 1)`, `9.0 u` at `f32`, flat in the scale; divided by the
    /// *result* instead, the same draws read `106 u` once the columns leave the `O(1)` range,
    /// which is a property of the question and not of this routine (`0048` decision 4).
    ///
    /// [`lminus`]: LieGroup::lminus
    #[inline]
    pub fn mul_inv(&self, other: &Self) -> Self {
        let q = self.q * other.q.conjugate();
        let r = SO3::from_quat_unchecked(q);
        Self {
            q,
            x: array::from_fn(|i| self.x[i] - r.act(other.x[i])),
        }
    }

    /// `other⁻¹ · self`, without forming `other⁻¹` (`0048`).
    ///
    /// `(q_b* q_a, R(q_b*)(x_{a,i} − x_{b,i}))`: **one** rotation per column, where
    /// `other.inverse() * *self` rotates `x_{b,i}` and `x_{a,i}` separately and subtracts after,
    /// so it does two and subtracts two rotated vectors instead of rotating one difference. This
    /// is the direction a relative transform `T_w_a⁻¹ · T_w_b` takes, and it is [`rminus`]'s
    /// first half — the **default** side (`0002`) — which is why that method's body is this one.
    ///
    /// Its own `NUMERICS.md` §14 twin, for [`mul_inv`](Self::mul_inv)'s reasons and with its
    /// scale-relative accuracy (`sen3_tests::group::sen3_inv_mul_matches_reference_n1`): measured
    /// worst `8.1 u` of `max(‖x_a‖, ‖x_b‖, 1)`, `8.0 u` at `f32`. Subtracting *before* the rotation is also the better-conditioned order — the
    /// cancellation happens in the inputs, where it is exact for nearby frames, rather than
    /// between two separately rounded rotations.
    ///
    /// [`rminus`]: LieGroup::rminus
    #[inline]
    pub fn inv_mul(&self, other: &Self) -> Self {
        // One conjugate, read twice: `SO3::inverse` *is* `Quat::conjugate`, so taking it through
        // `other.rotation().inverse()` as well would negate the same three components again.
        let qi = other.q.conjugate();
        let rinv = SO3::from_quat_unchecked(qi);
        Self {
            q: qi * self.q,
            x: array::from_fn(|i| rinv.act(self.x[i] - other.x[i])),
        }
    }

    /// The rotation part `R`.
    #[inline]
    pub fn rotation(&self) -> SO3<S> {
        SO3::from_quat_unchecked(self.q)
    }

    /// The element of a rotation and its columns, exactly: a move of `3 + 4` numbers, no
    /// transcendental and no round trip through `Exp`/`Log` (`0042`, which is `0029` one
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
    /// the only way to read an `N > 2` element without a `Log` (`0042`).
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
            true => (-mul_hat(&r, Vector(p.0)), r),
            false => (-hat(Vector(self.act(p).0)), Mat3::identity()),
        };
        let cols = array::from_fn(|c| match c < 3 {
            true => left.col(c),
            false => right.col(c - 3),
        });
        (Matrix::from_cols(cols), r)
    }

    /// `X⁻¹ p = Rᵀ (p − t)` (`NUMERICS.md` §5.5): the subtraction first, then
    /// [`SO3::inv_act`] (`0103`). `X.inverse() * p` rotates `p` and `t` apart and cancels two
    /// vectors of size `‖t‖` down to `‖X⁻¹ p‖`; here the difference rounds once relative to
    /// itself and the rotation sees its scale. It is
    /// [`inv_mul`](Self::inv_mul)'s column arithmetic, to the bit.
    #[inline]
    pub fn inv_act(&self, p: Point3<S>) -> Point3<S> {
        Point(self.rotation().inv_act(Vector(p.0) - self.x[0]).0)
    }

    /// `(∂(X⁻¹ p)/∂X, ∂(X⁻¹ p)/∂p)` in the convention of `Sd` (`NUMERICS.md` §2.4): rotation-first
    /// `[[X⁻¹ p]_×, −I]` on the right, `[Rᵀ [p]_×, −Rᵀ]` on the left, and `Rᵀ` (`0103`).
    ///
    /// The rotation block is [`SO3::inv_act_jacobians`]' at `p − t` (right) or `p` (left): the
    /// right one is `hat(self.inv_act(p))` to the bit, the left one `Rᵀ [p]_×` with no `t` to
    /// cancel. Dense for [`act_jacobians`](Self::act_jacobians)' reason (`0005`).
    #[inline]
    pub fn inv_act_jacobians<Sd: Side>(&self, p: Point3<S>) -> (Matrix<S, 3, 6>, Mat3<S>) {
        let v = match Sd::IS_RIGHT {
            true => Vector(p.0) - self.x[0],
            false => Vector(p.0),
        };
        let (rot, rt) = self.rotation().inv_act_jacobians::<Sd>(v);
        // `−I` written, not negated: `-Mat3::identity()` would put `−0` off its diagonal.
        let trans = match Sd::IS_RIGHT {
            true => Mat3::from_cols(array::from_fn(|c| {
                Vector(array::from_fn(|r| match r == c {
                    true => -S::one(),
                    false => S::zero(),
                }))
            })),
            false => -rt,
        };
        let cols = array::from_fn(|c| match c < 3 {
            true => rot.col(c),
            false => trans.col(c - 3),
        });
        (Matrix::from_cols(cols), rt)
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

/// `J_r(τ)` (`LEFT = false`) or `J_l(τ)` (`LEFT = true`), each column in SE.8 remark (ii)'s five
/// words (`NUMERICS.md` §5.3, `0097`, draft): `Q(ρ, φ) = S(ρ) + [aρ − (φ·ρ)(b − 2d)φ]^`, with the
/// symmetric part `S(ρ) = b(φ^ρ^ + ρ^φ^) − 2e(φ·ρ)φ^φ^` and `a = ½ − zd` SO(3)'s `(1 − cos θ)/θ²`
/// (CO.4(a)). `J_r`'s column `Q(−ρ, −φ)` has the same `S` and the skew part negated. The diagonal
/// is `SO3::jr(φ)` from the same coefficients, and its transpose is `SO3::jl(φ)` to the bit.
#[inline]
fn forward<S: Real, const N: usize, const LEFT: bool>(tau: &SEn3Tangent<S, N>) -> SEn3Jac<S, N> {
    let phi = tau.phi;
    let w2 = HatSq::of(phi);
    let ([a, b, t, e], sinc, far) = q_five_coeffs(norm_sq(phi));
    let right = w2.poly_far(S::one(), -a, b, sinc, far);
    let k = Cols::of(phi, b, -(e + e));
    SEn3Jac {
        diag: if LEFT { right.transpose() } else { right },
        col: tau.rho.map(|r| {
            let (sym, p) = k.sym(r);
            let h = r.scale(a) - phi.scale(p * t);
            sym.plus_hat(if LEFT { h } else { -h })
        }),
    }
}

/// `J_r⁻¹(τ)` (`LEFT = false`) or `J_l⁻¹(τ)` (`LEFT = true`) in `NUMERICS.md` §5.4's closed form.
///
/// The provided `jl_inv`, `J_r⁻¹(−τ)`, is this with `LEFT` set to the bit: the diagonal is
/// `SO3::jr_inv(−φ)`, which is `SO3::jr_inv(φ)ᵗ` exactly (`so3_tests::jl_is_jr_transposed_to_the_bit`),
/// and negating `φ` and `ρ` changes no product in [`InvCols::sym`] and only the sign of `½ρ`.
#[inline]
fn inverse<S: Real, const N: usize, const LEFT: bool>(tau: &SEn3Tangent<S, N>) -> SEn3Jac<S, N> {
    let (c, dc) = jr_inv_coeff_dz(norm_sq(tau.phi));
    let right = HatSq::of(tau.phi).poly(S::one(), S::lit(0.5), c);
    let k = Cols::of(tau.phi, c, dc + dc);
    SEn3Jac {
        diag: if LEFT { right.transpose() } else { right },
        col: tau.rho.map(|r| {
            let h = r.scale(S::lit(0.5));
            k.sym(r).0.plus_hat(if LEFT { -h } else { h })
        }),
    }
}

/// `(J_r⁻¹(τ), J_l⁻¹(τ))`, bit-identical to [`LieGroup::jr_inv`] and [`LieGroup::jl_inv`] called
/// separately, from one coefficient group and one symmetric part per column: the two sides differ
/// only in the diagonal's transpose and the sign of `½ρ^`.
#[inline]
fn inverses<S: Real, const N: usize>(tau: &SEn3Tangent<S, N>) -> (SEn3Jac<S, N>, SEn3Jac<S, N>) {
    let (c, dc) = jr_inv_coeff_dz(norm_sq(tau.phi));
    let right = HatSq::of(tau.phi).poly(S::one(), S::lit(0.5), c);
    let k = Cols::of(tau.phi, c, dc + dc);
    let mut jr = SEn3Jac {
        diag: right,
        col: [zero3(); N],
    };
    let mut jl = SEn3Jac {
        diag: right.transpose(),
        col: [zero3(); N],
    };
    for (i, &r) in tau.rho.iter().enumerate() {
        let (sym, h) = (k.sym(r).0, r.scale(S::lit(0.5)));
        jr.col[i] = sym.plus_hat(h);
        jl.col[i] = sym.plus_hat(-h);
    }
    (jr, jl)
}

/// The symmetric part `S(ρ) = c(φ^ρ^ + ρ^φ^) + c₂(φ·ρ)φ^φ^` of a column of `J(τ)` or `J⁻¹(τ)`:
/// `J⁻¹`'s with `c₂ = 2c′`, `c′ = dc/dθ²` (`NUMERICS.md` §5.4, `docs/maths/se3.md` SE.17), `J`'s
/// with `(b, −2e)` (§5.3, SE.8 remark (ii)). What depends on `φ` alone is formed once.
struct Cols<S> {
    phi: Vec3<S>,
    c: S,
    c2: S,
    /// `−(φ^φ^)ᵢᵢ`: `y² + z²`, `x² + z²`, `x² + y²`.
    sq: [S; 3],
    /// `(φ^φ^)ᵢⱼ`: `xy`, `xz`, `yz`.
    off: [S; 3],
}

impl<S: Real> Cols<S> {
    #[inline]
    fn of(phi: Vec3<S>, c: S, c2: S) -> Self {
        let [x, y, z] = phi.0;
        Self {
            phi,
            c,
            c2,
            sq: [z * z + y * y, z * z + x * x, y * y + x * x],
            off: [x * y, x * z, y * z],
        }
    }

    /// `(S(ρ), φ·ρ)`. `φ^ρ^ + ρ^φ^ = φρᵀ + ρφᵀ − 2(φ·ρ)I`, so each diagonal entry of `S` is
    /// `−(2c(φⱼρⱼ + φₖρₖ) + c₂(φ·ρ)(φⱼ² + φₖ²))`: the `(φ·ρ)` and `θ²` that would cancel against
    /// `φᵢρᵢ` and `φᵢ²` are never formed, as `HatSq` forms `W²`'s diagonal.
    #[inline]
    fn sym(&self, rho: Vec3<S>) -> (Sym<S>, S) {
        let [x, y, z] = self.phi.0;
        let [r0, r1, r2] = rho.0;
        let p = self.phi.dot(rho);
        let q = self.c2 * p;
        let cc = self.c + self.c;
        let (px, py, pz) = (x * r0, y * r1, z * r2);
        let [s0, s1, s2] = self.sq;
        let [o01, o02, o12] = self.off;
        let sym = Sym {
            diag: [
                -(cc * (py + pz) + q * s0),
                -(cc * (px + pz) + q * s1),
                -(cc * (px + py) + q * s2),
            ],
            off: [
                self.c * (x * r1 + r0 * y) + q * o01,
                self.c * (x * r2 + r0 * z) + q * o02,
                self.c * (y * r2 + r1 * z) + q * o12,
            ],
        };
        (sym, p)
    }
}

/// A symmetric `3 x 3`: its diagonal and its `(0, 1)`, `(0, 2)`, `(1, 2)` entries.
#[derive(Clone, Copy)]
struct Sym<S> {
    diag: [S; 3],
    off: [S; 3],
}

impl<S: Real> Sym<S> {
    /// `self + h^`, entry by entry: the hat's zeros are not added.
    #[inline]
    fn plus_hat(&self, h: Vec3<S>) -> Mat3<S> {
        let [h0, h1, h2] = h.0;
        let [d0, d1, d2] = self.diag;
        let [s01, s02, s12] = self.off;
        Matrix::from_rows([
            Vector([d0, s01 - h2, s02 + h1]),
            Vector([s01 + h2, d1, s12 - h0]),
            Vector([s02 - h1, s12 + h0, d2]),
        ])
    }
}

/// `PHASE4.md` §1.2's screw geodesic, in `SO3::geodesic`'s two arms (`0054`) and in the **world
/// frame**: conjugating `Δ = X₀⁻¹X₁` by `R₀` makes its translation `x₁ − x₀` and its axis
/// `v_w = vec(q₁ q₀*)` (GE.15), so no column is rotated in or out. The rotation is
/// [`geodesic_short`]'s or [`geodesic_long`]'s, `SO3::geodesic`'s to the bit, and the hot path is
/// one branch: the short arm forms rotation and translation together.
#[inline]
fn screw_geodesic<S: Real, const N: usize>(x0: &SEn3<S, N>, x1: &SEn3<S, N>, t: S) -> SEn3<S, N> {
    let rel = geodesic_rel(x0.q, x1.q, t);
    let m = x1.q * x0.q.conjugate();
    let v = Vector([rel.flip * m.x, rel.flip * m.y, rel.flip * m.z]);
    let dx: [Vec3<S>; N] = array::from_fn(|i| x1.x[i] - x0.x[i]);
    let (q, step) = S::branch(
        rel.fast,
        || {
            let (rot, r) = geodesic_short(x0.q, &rel, t);
            (rot, screw_short(v, r, t, &dx))
        },
        || {
            let (rot, s1, cos_t, nv) = geodesic_long::<S, true>(x0.q, x1.q, &rel, t);
            // Exactly `t = 1` below the switch is the one way here with `short` set; its
            // translation is the short arm's, at `log_ratio_short`'s safe argument by `short`
            // itself.
            let step = S::branch(
                rel.short,
                || screw_short(v, log_ratio_short(rel.n2, rel.w, rel.short), t, &dx),
                || screw_long(v, &rel, nv, s1, cos_t, t, &dx),
            );
            (rot, step)
        },
    );
    SEn3 {
        q,
        x: array::from_fn(|i| x0.x[i] + step[i]),
    }
}

/// The screw twin's translation below `r`'s second switch: the definition,
/// `J_l(tφ) t J_l⁻¹(φ) (x₁ − x₀)` with `φ = r v_w`, by two cross products each way as `log` and
/// `exp` write them. Every coefficient is the catalogue's on its series arm: no transcendental,
/// and through `Dual` the provided body's accuracy.
#[inline]
fn screw_short<S: Real, const N: usize>(v: Vec3<S>, r: S, t: S, dx: &[Vec3<S>; N]) -> [Vec3<S>; N] {
    let half = S::lit(0.5);
    let phi = v.scale(r);
    let c = jr_inv_coeff(norm_sq(phi));
    let phit = phi.scale(t);
    let (a, b) = jr_coeffs(norm_sq(phit));
    dx.map(|d| {
        let c1 = phi.cross(d);
        let rho = ((d - c1.scale(half)) + phi.cross(c1).scale(c)).scale(t);
        let c2 = phit.cross(rho);
        (rho + c2.scale(a)) + phit.cross(c2).scale(b)
    })
}

/// The screw twin's translation above `r`'s second switch: GE.12's dual part,
/// `q_d = ½(0, x₁ − x₀) ⊗ (w, v_w)`, `ϰ = q_d,w/‖v‖²`, `m̄ = q_d,v + ϰ w v_w`,
/// `q_dᵗ = (t ϖ_t q_d,w, ϖ_t m̄ − t ϰ cos tα v_w)`, `2 vec(q_dᵗ (cos tα, ϖ_t v_w)*)`, with
/// `ϖ_t = sin(tα)/‖v‖` -- not the blend's `sin(tα)/sin α`, which errs more on every regime
/// `0054` measured. `ϰ` carries GE.13(c)'s `10² α⁻¹ u` derivative loss, which the switch bounds,
/// and `‖v‖²` is never `0` where this arm is selected; the `select`s are `0003`'s safe arguments
/// for a lane that is not.
#[inline]
fn screw_long<S: Real, const N: usize>(
    v: Vec3<S>,
    rel: &crate::so3::GeodesicRel<S>,
    nv: S,
    s1: S,
    cos_t: S,
    t: S,
    dx: &[Vec3<S>; N],
) -> [Vec3<S>; N] {
    let (one, half, two) = (S::one(), S::lit(0.5), S::lit(2.0));
    let point = nv.le(S::zero());
    let varpi = S::select(point, t, s1 / S::select(point, one, nv));
    let n2 = S::select(point, one, rel.n2);
    let (w, [vx, vy, vz]) = (rel.w, v.0);
    let power_conj = Quat {
        w: cos_t,
        x: -(varpi * vx),
        y: -(varpi * vy),
        z: -(varpi * vz),
    };
    dx.map(|d| {
        // `½(0, d) ⊗ (w, v)`, the Hamilton product without its zero lanes: the same bits as
        // `Quat`'s `Mul` but for the sign of an exact zero, and `screw_pow`'s grouping.
        let [dx, dy, dz] = d.0;
        let dw = half * (((-(dx * vx)) - dy * vy) - dz * vz);
        let dv = Vector([
            (dx * w + dy * vz) - dz * vy,
            ((-(dx * vz)) + dy * w) + dz * vx,
            (dx * vy - dy * vx) + dz * w,
        ])
        .scale(half);
        let kappa = dw / n2;
        let mbar = dv + v.scale(kappa * w);
        let [ex, ey, ez] = (mbar.scale(varpi) - v.scale((t * kappa) * cos_t)).0;
        let xt = Quat {
            w: (t * varpi) * dw,
            x: ex,
            y: ey,
            z: ez,
        } * power_conj;
        Vector([two * xt.x, two * xt.y, two * xt.z])
    })
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
    /// what `tf_tree_math::exp_se3` does below `b`'s switch (`0010`, `0041` draft).
    ///
    /// Where `b` is on its closed arm, `ρ + b (φ × (φ × ρ))` cancels to `(sin θ/θ) ρ` near `π`, so
    /// the column is `(sin θ/θ) ρ + a (φ × ρ) + b (φ·ρ) φ` instead, with the `sin θ/θ` that arm
    /// formed `b` from: `J_r`'s diagonal of `0089`, applied (`0095`, draft). One branch on `b`'s
    /// own mask, which `jr_coeffs_sinc` has already taken; the quaternion is formed first, which
    /// measured 4–7 ns faster at `f64` than after it.
    #[inline]
    fn exp(tau: &SEn3Tangent<S, N>) -> Self {
        let phi = tau.phi;
        let (a, b, sinc, far) = jr_coeffs_sinc(norm_sq(phi));
        let q = SO3::exp(&SO3Tangent { phi }).quat();
        let rho = tau.rho;
        let x = S::branch(
            far,
            move || {
                rho.map(|p| (p.scale(sinc) + phi.cross(p).scale(a)) + phi.scale(b * phi.dot(p)))
            },
            move || {
                rho.map(|p| {
                    let c1 = phi.cross(p);
                    (p + c1.scale(a)) + phi.cross(c1).scale(b)
                })
            },
        );
        Self { q, x }
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

    /// `Log(other⁻¹ · self)` through [`inv_mul`](Self::inv_mul) (`0048` decision 3).
    ///
    /// The provided body is `(base.inverse() * *self).log()`, which rotates two columns per `i`
    /// where `inv_mul` rotates one: `3N` `act`s become `N`, and the subtraction moves ahead of
    /// the rotation, where it is exact for nearby frames. Same map, fewer roundings.
    #[inline]
    fn rminus(&self, base: &Self) -> SEn3Tangent<S, N> {
        self.inv_mul(base).log()
    }

    /// `Log(self · other⁻¹)` through [`mul_inv`](Self::mul_inv), for
    /// [`rminus`](Self::rminus)'s reasons.
    #[inline]
    fn lminus(&self, base: &Self) -> SEn3Tangent<S, N> {
        self.mul_inv(base).log()
    }

    /// SE(3): `PHASE4.md` §1.2's dual-quaternion power, `screw_geodesic` (`0054`). Other `N`
    /// keep the provided body, since no `se23_geodesic` stratum verifies the twin there (`0006`).
    #[inline]
    fn geodesic(x0: &Self, x1: &Self, t: S) -> Self {
        if N == 1 {
            screw_geodesic(x0, x1, t)
        } else {
            crate::reference::geodesic(x0, x1, t)
        }
    }

    /// `Ad_X = R + ε [x_i]_× R` (`NUMERICS.md` §5.2).
    #[inline]
    fn adjoint(&self) -> SEn3Jac<S, N> {
        let r = self.rotation().to_matrix();
        SEn3Jac {
            diag: r,
            col: self.x.map(|v| hat_mul(v, &r)),
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

    /// `J_r(τ) = J_r(φ) + ε Q(−ρ_i, −φ)` (`NUMERICS.md` §5.3), each column in five words: `forward`.
    #[inline]
    fn jr(tau: &SEn3Tangent<S, N>) -> SEn3Jac<S, N> {
        forward::<S, N, false>(tau)
    }

    /// `J_l(τ) = J_l(φ) + ε Q(ρ_i, φ)`: the provided `J_r(−τ)` to the bit, without negating the
    /// tangent.
    #[inline]
    fn jl(tau: &SEn3Tangent<S, N>) -> SEn3Jac<S, N> {
        forward::<S, N, true>(tau)
    }

    /// `J_r⁻¹(τ) = J_r⁻¹(φ) + ε(½ρᵢ^ + c(φ^ρᵢ^ + ρᵢ^φ^) + 2c′(φ·ρᵢ)φ^φ^)` (`NUMERICS.md` §5.4,
    /// `0093`, draft).
    ///
    /// The dual-matrix inverse of [`jr`](LieGroup::jr), `J_r⁻¹(φ) − ε J_r⁻¹(φ) Q(−ρᵢ, −φ) J_r⁻¹(φ)`,
    /// reduced to four words (`docs/maths/se3.md` SE.17): no `Q`, no `3 x 3` product, and `c′` the
    /// `Dual` lane of `c`'s own evaluation. Against the sandwich, 4.2–6.8 `u` → 1.4 `u` at worst on
    /// the corpus, and 3–5× faster. Like the sandwich it divides by nothing, so §12's finite value
    /// holds; [`Jac::inverse`] remains the inverse of a `SEn3Jac` that is nobody's `J_r` (`0005`).
    ///
    /// # Domain
    ///
    /// `θ < 2π`, [`SO3::jr_inv`]'s, asserted nowhere for its reason.
    #[inline]
    fn jr_inv(tau: &SEn3Tangent<S, N>) -> SEn3Jac<S, N> {
        inverse::<S, N, false>(tau)
    }

    /// `J_l⁻¹(τ)`: the provided `J_r⁻¹(−τ)` to the bit, without negating the tangent.
    #[inline]
    fn jl_inv(tau: &SEn3Tangent<S, N>) -> SEn3Jac<S, N> {
        inverse::<S, N, true>(tau)
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
    /// do not: `J_l⁻¹(τ)`'s carries `−½ρ^` and `J_r⁻¹(τ)`'s `+½ρ^`, and the transpose of
    /// a `SEn3Jac` is block *upper* triangular, so it is not a `SEn3Jac` to return.
    ///
    /// Fused through this module's `inverses`: one coefficient group and one symmetric part per
    /// column serve both sides. `the_fused_inverses_are_the_separate_ones_to_the_bit` holds it to
    /// the unfused pair, and
    /// `laws::jacobian_rows` holds this row against separate `jr_inv`/`jl_inv` calls at a bound of
    /// exactly `0`.
    #[inline]
    fn rminus_jacobians(&self, base: &Self) -> (SEn3Jac<S, N>, SEn3Jac<S, N>) {
        let (jr, jl) = inverses(&self.rminus(base));
        (jr, jl.neg())
    }

    /// `(J_l⁻¹(τ), −J_r⁻¹(τ))` at `τ = self ⊖_L base` (`NUMERICS.md` §2.3); see
    /// [`rminus_jacobians`](LieGroup::rminus_jacobians) for the two inversions.
    #[inline]
    fn lminus_jacobians(&self, base: &Self) -> (SEn3Jac<S, N>, SEn3Jac<S, N>) {
        let (jr, jl) = inverses(&self.lminus(base));
        (jl, jr.neg())
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
