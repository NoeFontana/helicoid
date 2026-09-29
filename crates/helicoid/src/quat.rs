//! The quaternion `Quat` (`docs/PHASE3.md` §4): Hamilton, `w` first, active (`0002`).

use core::ops::Mul;
use helicoid_linalg::{Blend, Mask, Mat3, Precision, Real, Vector};

/// A quaternion `w + x i + y j + z k`, Hamilton (`i² = j² = k² = ijk = -1`), stored `w` first.
///
/// Active: a unit `q` rotates `v` to `q v q*` ([`to_matrix`](Quat::to_matrix), `NUMERICS.md` §1).
/// `Quat` is a value with no invariant: the unit-norm contract belongs to the group type built on
/// it, and each constructor states what it does about it. Composition is the Hamilton product,
/// spelled `p * q` (`docs/API.md` R1), and never normalizes (`NUMERICS.md` §3.6); there is no
/// `Add` or `Sub`.
///
/// Other component orders enter and leave only through named converters (R3): [`from_xyzw`],
/// [`to_xyzw`] and [`from_jpl`].
///
/// [`from_xyzw`]: Quat::from_xyzw
/// [`to_xyzw`]: Quat::to_xyzw
/// [`from_jpl`]: Quat::from_jpl
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct Quat<S> {
    /// The scalar part.
    pub w: S,
    /// The coefficient of `i`.
    pub x: S,
    /// The coefficient of `j`.
    pub y: S,
    /// The coefficient of `k`.
    pub z: S,
}

/// `2^-40` (`f64`) or `2^-16` (`f32`): the bound of `NUMERICS.md` §3.6 on `| ‖q‖² - 1 |`.
///
/// Both are exactly representable, so `lit` is exact; the division is exact too.
fn unit_tolerance<S: Real>() -> S {
    match S::PRECISION {
        Precision::F64 => S::lit(1.0 / 1_099_511_627_776.0),
        Precision::F32 => S::lit(1.0 / 65_536.0),
    }
}

impl<S: Real> Quat<S> {
    /// The quaternion `(w, x, y, z)`, taken as it is.
    ///
    /// # Domain
    ///
    /// `| ‖q‖² - 1 | <= 2^-40` (`f64`, `Dual<f64, N>`) or `2^-16` (`f32`) (`NUMERICS.md` §3.6,
    /// §12), checked by `debug_assert!` (a NaN fails it); the bound is chosen by
    /// [`Real::PRECISION`]. Release builds check nothing and return the four values, so a
    /// quaternion outside the domain is garbage in, garbage out. One [`renormalize`] step repairs
    /// a quaternion inside it to rounding level (`docs/maths/so3.md` SO.15).
    ///
    /// [`renormalize`]: Quat::renormalize
    #[inline]
    pub fn from_wxyz_unchecked(w: S, x: S, y: S, z: S) -> Self {
        let q = Self { w, x, y, z };
        debug_assert!(
            (q.norm_sq() - S::one())
                .abs()
                .le(unit_tolerance::<S>())
                .all(),
            "Quat::from_wxyz_unchecked: | |q|^2 - 1 | exceeds 2^-40 (f64) or 2^-16 (f32)"
        );
        q
    }

    /// The quaternion `(w, x, y, z)` after one [`renormalize`](Quat::renormalize) step.
    ///
    /// This is the first-order Newton step, not a normalization: it is exact to
    /// `O((‖q‖² - 1)²)` and no more, and it does not compute a `sqrt`. Its accuracy domain is
    /// `| ‖q‖² - 1 | <= 2^-26.29` (`f64`) or `2^-11.79` (`f32`) (`docs/maths/so3.md` SO.15);
    /// see [`renormalize`](Quat::renormalize).
    ///
    /// # Domain
    ///
    /// None is asserted (`NUMERICS.md` §3.6 and §12 state none); a release and a debug build
    /// give the same result for every input.
    #[inline]
    pub fn from_wxyz_normalized(w: S, x: S, y: S, z: S) -> Self {
        let mut q = Self { w, x, y, z };
        q.renormalize();
        q
    }

    /// The quaternion of a `[x, y, z, w]` array (Eigen, nalgebra's `coords`, ROS): a permutation,
    /// nothing else. The norm is not checked, as for a struct literal.
    #[inline]
    pub fn from_xyzw(a: [S; 4]) -> Self {
        Self {
            w: a[3],
            x: a[0],
            y: a[1],
            z: a[2],
        }
    }

    /// The components as `[x, y, z, w]`, the inverse of [`from_xyzw`](Quat::from_xyzw).
    #[inline]
    pub fn to_xyzw(&self) -> [S; 4] {
        [self.x, self.y, self.z, self.w]
    }

    /// The Hamilton quaternion `(w, -x, -y, -z)` of a JPL quaternion given as `[x, y, z, w]`.
    ///
    /// JPL quaternions (Trawny and Roumeliotis 2005) multiply in the opposite order and their
    /// rotation matrix, for a unit `q`, is `C(q) = (2w² - 1) I - 2w [u]× + 2 u uᵀ`, the transpose
    /// of the Hamilton `R` of the same four numbers. The conjugate has exactly that matrix:
    /// [`to_matrix`] of the result is `C(q)` for a unit `q` (Sommer et al. 2018) and `Rᵀ` for
    /// every `q`. Which frame the matrix maps from and to stays the caller's: this is a change of
    /// representation, and no transposition is applied on top.
    ///
    /// [`to_matrix`]: Quat::to_matrix
    #[inline]
    pub fn from_jpl(a: [S; 4]) -> Self {
        Self {
            w: a[3],
            x: -a[0],
            y: -a[1],
            z: -a[2],
        }
    }

    /// The conjugate `(w, -x, -y, -z)`; the inverse of a unit quaternion.
    #[inline]
    pub fn conjugate(&self) -> Self {
        Self {
            w: self.w,
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }

    /// `‖q‖² = w² + x² + y² + z²`, summed in that order. No `sqrt`.
    #[inline]
    pub fn norm_sq(&self) -> S {
        self.w * self.w + self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// The matrix `R(q) = (w² - ‖u‖²) I + 2 u uᵀ + 2w [u]×` of `NUMERICS.md` §1, `u = (x, y, z)`.
    ///
    /// For a unit `q` it is the rotation `v ↦ q v q*`, orthogonal with determinant `+1`. For any
    /// other `q` it is `‖q‖² R(q / ‖q‖)` (`docs/maths/so3.md` SO.2), a scaled rotation.
    ///
    /// # Domain
    ///
    /// None is asserted. The rotation reading needs `q` unit to the tolerance of
    /// [`from_wxyz_unchecked`](Quat::from_wxyz_unchecked): the matrix errs by
    /// `| ‖q‖² - 1 | ‖v‖` when applied to `v`.
    pub fn to_matrix(&self) -> Mat3<S> {
        let (w, x, y, z) = (self.w, self.x, self.y, self.z);
        let two = S::lit(2.0);
        let (tw, tx, ty, tz) = (two * w, two * x, two * y, two * z);
        let d = w * w - (x * x + y * y + z * z);
        // Column-major: `col(c)[r]`; `[u]×` is `[[0, -z, y], [z, 0, -x], [-y, x, 0]]`.
        Mat3::from_cols([
            Vector([d + tx * x, tx * y + tw * z, tx * z - tw * y]),
            Vector([tx * y - tw * z, d + ty * y, ty * z + tw * x]),
            Vector([tx * z + tw * y, ty * z - tw * x, d + tz * z]),
        ])
    }

    /// One first-order Newton step, `q ← q (3 - ‖q‖²)/2` (`NUMERICS.md` §3.6), in place. It is
    /// the only `&mut self` numeric method (`docs/API.md` R2).
    ///
    /// With `η = ‖q‖² - 1` the step gives `‖q'‖² - 1 = -¾η² + ¼η³` exactly and keeps the
    /// direction of `q` for `η < 2` (`docs/maths/so3.md` SO.14). The truncation term `¾η²` is at
    /// most `u` for `|η| <= 2^-26.29` (`f64`) or `2^-11.79` (`f32`); the computed `|η'|` is a few
    /// `u` (measured at most `5u` for `f64` and `3.73u` for `f32`, from the edge of the domain of
    /// [`from_wxyz_unchecked`](Quat::from_wxyz_unchecked)). From a larger `η` it is only a better
    /// guess, not a normalization: it returns zero at `η = 2` and reverses `q` beyond it. Divide
    /// by `sqrt(‖q‖²)` for a full normalization.
    ///
    /// # Domain
    ///
    /// None is asserted (`NUMERICS.md` §3.6 and §12 state none), so the step is defined for every
    /// input, as `docs/maths/so3.md` SO.14 reads it; only its accuracy has a domain, above.
    #[inline]
    pub fn renormalize(&mut self) {
        let k = (S::lit(3.0) - self.norm_sq()) * S::lit(0.5);
        *self = Self {
            w: self.w * k,
            x: self.x * k,
            y: self.y * k,
            z: self.z * k,
        };
    }

    // The signed sums live here so that clippy's `suspicious_arithmetic_impl` does not read the
    // `-` in `Mul::mul` as a typo (as in `Rn`). Each component is a left-to-right sum of four
    // products (`NUMERICS.md` §1: `(w,u)(w',u') = (ww' - u·u', wu' + w'u + u×u')`).
    #[inline]
    fn hamilton(self, o: Self) -> Self {
        Self {
            w: self.w * o.w - self.x * o.x - self.y * o.y - self.z * o.z,
            x: self.w * o.x + self.x * o.w + self.y * o.z - self.z * o.y,
            y: self.w * o.y - self.x * o.z + self.y * o.w + self.z * o.x,
            z: self.w * o.z + self.x * o.y - self.y * o.x + self.z * o.w,
        }
    }
}

/// The Hamilton product: `a * b` is `T_a_x · T_x_b` on unit quaternions (`0002`).
///
/// It does not normalize; the norm drifts by at most `~16u` per product
/// (`docs/maths/so3.md` SO.15).
impl<S: Real> Mul for Quat<S> {
    type Output = Self;
    #[inline]
    fn mul(self, o: Self) -> Self {
        self.hamilton(o)
    }
}

impl<S: Real> Blend<S> for Quat<S> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self {
            w: S::select(m, t.w, f.w),
            x: S::select(m, t.x, f.x),
            y: S::select(m, t.y, f.y),
            z: S::select(m, t.z, f.z),
        }
    }
}
