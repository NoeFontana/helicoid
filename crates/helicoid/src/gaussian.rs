//! A concentrated Gaussian on a Lie group, its side in the type (`docs/PHASE5.md` §5, `0065`).
//!
//! `Gaussian<S, G, Right, D>` is the law of `μ Exp(ξ)` and `Gaussian<S, G, Left, D>` that of
//! `Exp(ξ) μ`, `ξ ~ N(0, Σ)` in rotation-first coordinates (`docs/maths/gamma-gaussian.md` GG.7).
//! The side is a type parameter because a covariance read on the wrong side is the most expensive
//! of GG.14's seven misreadings (V1, `NEES = 3.2·10³` at a 197 m pose) and nothing in a raw matrix
//! shows it.

use crate::{Jac, Left, LieGroup, Right, Side, Tangent};
use core::array;
use core::marker::PhantomData;
use helicoid_linalg::{chol, solve_lower, Matrix, Real, Vector};

/// The law of `mean ⊕_Sd ξ`, `ξ ~ N(0, cov)` (GG.7): `mean` is a reference point with `E ξ = 0`,
/// not a Fréchet mean, and the definition is meant for a concentrated law (`3σ_φ ≪ π`).
///
/// **`cov` is exactly symmetric**: [`new`](Self::new), the only way in, keeps the lower triangle
/// and copies it over the upper one, and `cov` is read through [`cov`](Self::cov) (`0066`). A
/// `J Σ Jᵀ` formed in floating point is not symmetric (up to `1.2u` of its largest entry, GG.11),
/// and a consumer that factors or eigen-decomposes it would otherwise see two matrices. The copy
/// is `D(D−1)/2` moves and no arithmetic, paid once per value rather than on every read.
///
/// `D` is `G::DOF`, a parameter only because stable Rust cannot size an array from an associated
/// const; [`new`](Self::new) asserts it at monomorphization:
///
/// ```compile_fail
/// use helicoid::{Gaussian, LieGroup, Right, SE3};
/// use helicoid_linalg::Matrix;
/// let _ = Gaussian::<f64, SE3<f64>, Right, 5>::new(SE3::identity(), Matrix::identity());
/// ```
///
/// Positive control:
///
/// ```
/// use helicoid::{Gaussian, LieGroup, Right, SE3};
/// use helicoid_linalg::Matrix;
/// let g = Gaussian::<f64, SE3<f64>, Right, 6>::new(SE3::identity(), Matrix::identity());
/// let (d2, ok) = g.mahalanobis_sq(&SE3::identity());
/// assert!(ok && d2 == 0.0);
/// ```
#[derive(Clone, Copy, Debug)]
pub struct Gaussian<S: Real, G: LieGroup<S>, Sd: Side, const D: usize> {
    /// The reference point `μ`.
    pub mean: G,
    /// `Σ`, rotation-first, in `Sd`'s convention, exactly symmetric; private so that holds.
    cov: Matrix<S, D, D>,
    _side: PhantomData<Sd>,
}

/// `D == G::DOF == G::Tangent::DOF`: the first sizes `cov`, the second the dense residual.
const fn tie<S: Real, G: LieGroup<S>, const D: usize>() {
    assert!(D == <G as LieGroup<S>>::DOF && D == <G::Tangent as Tangent<S>>::DOF);
}

/// `m` with its upper triangle replaced by its lower one.
fn lower_mirrored<S: Real, const D: usize>(m: &Matrix<S, D, D>) -> Matrix<S, D, D> {
    Matrix::from_cols(array::from_fn(|c| {
        Vector(array::from_fn(|r| {
            if r >= c {
                m.get(r, c)
            } else {
                m.get(c, r)
            }
        }))
    }))
}

impl<S: Real, G: LieGroup<S>, Sd: Side, const D: usize> Gaussian<S, G, Sd, D> {
    /// The law of `mean ⊕_Sd ξ`, `ξ ~ N(0, cov)`: **only `cov`'s lower triangle is read**, as
    /// [`chol`] reads it, and it is copied over the upper one.
    ///
    /// Nothing is asserted about `cov`: positive definiteness is what
    /// [`mahalanobis_sq`](Self::mahalanobis_sq)'s mask reports, and a value function owes no other
    /// check (`API.md` R4).
    #[inline]
    pub fn new(mean: G, cov: Matrix<S, D, D>) -> Self {
        const { tie::<S, G, D>() };
        Self {
            mean,
            cov: lower_mirrored(&cov),
            _side: PhantomData,
        }
    }

    /// `Σ`, rotation-first, in `Sd`'s convention: exactly symmetric.
    #[inline]
    pub fn cov(&self) -> &Matrix<S, D, D> {
        &self.cov
    }

    /// `(w, mask)`: `Σ = L Lᵀ` factored once by [`chol`], for testing many points against this
    /// law (`0066`). [`Whitener::mahalanobis_sq`] is then `⊖` and a forward substitution per point,
    /// `D²/2` multiply-adds and `D` divisions where `mahalanobis_sq` adds `chol`'s `D³/6` and `D`
    /// square roots, and is [`mahalanobis_sq`](Self::mahalanobis_sq)'s value to the bit.
    ///
    /// The mask is `chol`'s, and what it says is `mahalanobis_sq`'s: **a clear mask means no `d²`
    /// from `w` means anything** (GG.11(c)).
    #[inline]
    pub fn whitener(&self) -> (Whitener<S, G, Sd, D>, S::Mask) {
        const { tie::<S, G, D>() };
        let (l, mask) = chol(&self.cov);
        (
            Whitener {
                mean: self.mean,
                l,
                _side: PhantomData,
            },
            mask,
        )
    }

    /// `J Σ Jᵀ` about `mean` (GG.9(a)): the first-order law of `F(X)` for `X` this law and
    /// `j = D^{Sd} F(μ)`, `mean = F(μ)`.
    ///
    /// **`j` must be the Jacobian on side `Sd`** — `D^R F` for `Right`, `D^L F` for `Left`. The
    /// type cannot check it: `G::Jac` is one type for `J_r`, `J_l` and `Ad`. A Jacobian of the
    /// other side converts the output covariance by `Ad_{F(μ)}` and the input one by `Ad_μ⁻¹`,
    /// silently (GG.9; V7 of GG.14 is `J_l` for `J_r`, `NEES = 2.5·10³`). The mean is right to
    /// `O(Σ)`, the covariance to relative `O(Σ)` (GG.9(b)).
    #[inline]
    pub fn propagate(&self, j: &G::Jac, mean: G) -> Self {
        Self::new(mean, j.sandwich(&self.cov))
    }

    /// `(d², mask)`, `d² = ‖L⁻¹ (x ⊖_Sd mean)‖²` with `Σ = L Lᵀ` from [`chol`] (GG.10): a
    /// forward substitution, `Σ⁻¹` never formed, the squares summed in index order.
    ///
    /// The mask is `chol`'s, returned unaltered. **A clear mask means `d²` means nothing** — not
    /// that it is large, and not even that it is finite (GG.11(c)); a gate rejects on a clear mask
    /// rather than testing `d²` for finiteness. A set mask certifies the correlation matrix
    /// positive definite to about `n²u`, and `d²` is then accurate to `c_n u/λ_min(H)` relative,
    /// `c_6 ≈ 120` (GG.11–GG.12), whatever the units.
    ///
    /// `d²` is side-invariant in exact arithmetic (GG.10(a)), not in floating point: the left form
    /// of an SE(3) covariance at a large translation is the ill-conditioned one
    /// (`λ_min(H_L) ≈ σ_ρ²/(2σ_φ²‖t‖²)`, GG.13(b)), so a `χ²` gate is best computed on the
    /// right side.
    ///
    /// # Domain
    ///
    /// `θ(mean⁻¹ x) < π` (`Left`: `θ(x mean⁻¹) < π`), where `⊖` has its cut and `d²` jumps
    /// (GG.10(d)). Nothing is asserted.
    #[inline]
    pub fn mahalanobis_sq(&self, x: &G) -> (S, S::Mask) {
        const { tie::<S, G, D>() };
        // Not `whitener()` then its method: building the `Whitener` copies `L`, which measured
        // 1.3x to 1.4x this call's time, and `⊖` before `chol`, as `0065` shipped it: after, it
        // measured 1.05x to 1.16x (`0066`). `Whitener::mahalanobis_sq` is the same two lines.
        let delta = residual::<S, G, Sd, D>(&self.mean, x);
        let (l, mask) = chol(&self.cov);
        (solve_lower(&l, delta).norm_sq(), mask)
    }
}

impl<S: Real, G: LieGroup<S>, const D: usize> Gaussian<S, G, Right, D> {
    /// The same law as a left Gaussian: `Σ_L = Ad_μ Σ_R Ad_μᵀ`, **exactly**, not to first order
    /// (GG.8(a)), so the only error is rounding — which at a large translation is not small on
    /// the translation block (GG.13(c)).
    #[inline]
    pub fn to_left(&self) -> Gaussian<S, G, Left, D> {
        Gaussian::new(self.mean, self.mean.adjoint().sandwich(&self.cov))
    }
}

impl<S: Real, G: LieGroup<S>, const D: usize> Gaussian<S, G, Left, D> {
    /// The same law as a right Gaussian: `Σ_R = Ad_{μ⁻¹} Σ_L Ad_{μ⁻¹}ᵀ` (GG.8(a)), with the adjoint
    /// of the computed inverse rather than an inverted `Jac`: `Ad` of a group element is formed
    /// from its entries, where `Jac::inverse` would round a second time.
    #[inline]
    pub fn to_right(&self) -> Gaussian<S, G, Right, D> {
        Gaussian::new(self.mean, self.mean.inverse().adjoint().sandwich(&self.cov))
    }
}

/// A [`Gaussian`] with `Σ = L Lᵀ` factored, from [`Gaussian::whitener`]: a `χ²` gate's many
/// points against one law pay for `chol` once (`0066`).
#[derive(Clone, Copy, Debug)]
pub struct Whitener<S: Real, G: LieGroup<S>, Sd: Side, const D: usize> {
    mean: G,
    l: Matrix<S, D, D>,
    _side: PhantomData<Sd>,
}

impl<S: Real, G: LieGroup<S>, Sd: Side, const D: usize> Whitener<S, G, Sd, D> {
    /// `d² = ‖L⁻¹ (x ⊖_Sd mean)‖²` (GG.10): [`Gaussian::mahalanobis_sq`]'s value to the bit,
    /// without its mask, which [`Gaussian::whitener`] returned once, and with its accuracy.
    ///
    /// # Domain
    ///
    /// [`Gaussian::mahalanobis_sq`]'s. Nothing is asserted.
    #[inline]
    pub fn mahalanobis_sq(&self, x: &G) -> S {
        solve_lower(&self.l, residual::<S, G, Sd, D>(&self.mean, x)).norm_sq()
    }
}

/// `x ⊖_Sd mean`, dense and rotation-first.
#[inline]
fn residual<S: Real, G: LieGroup<S>, Sd: Side, const D: usize>(mean: &G, x: &G) -> Vector<S, D> {
    let mut delta = [S::zero(); D];
    Sd::minus(x, mean).write_dense(&mut delta);
    Vector(delta)
}
