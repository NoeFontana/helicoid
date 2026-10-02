//! The SE_N(3) tangent `SEn3Tangent` and the twist converters (`docs/PHASE3.md` §5).

use crate::traits::Tangent;
use core::array;
use core::iter::once;
use helicoid_linalg::{Blend, Real, Vec3, Vector};

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
