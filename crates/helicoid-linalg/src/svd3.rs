//! Singular value decomposition of a 3x3 matrix (`docs/PHASE2.md` §6, decision `0024`).
//!
//! A cyclic one-sided Hestenes Jacobi: right rotations orthogonalize the columns of `B = A V`, so
//! the singular values are the column norms and the left vectors the normalized columns. `A^T A`
//! is never formed, which is the whole reason `0024` amended §6 away from McAdams. The sweep count
//! is fixed, every decision is a mask, and every `sqrt` and division sees a safe argument, so all
//! lanes stay finite and a `Dual` result differentiates the code that ran. The gauge
//! (`det U = det V = +1`, the reflection in `sigma_3`) is built rather than repaired: no
//! determinant is evaluated anywhere.

use crate::eig3::{any_orthogonal, unit_or};
use crate::matrix::Mat3;
use crate::real::{Blend, Real};
use crate::vector::{Vec3, Vector};

/// Sweeps of the cycle `(0,1), (0,2), (1,2)`, fixed at every precision (`0024` decision 3).
///
/// Four. Over 20507 matrices in 17 strata (`f64`), the worst
/// `max_{j<k} |b_j . b_k| / sigma_1^2` is `4.0e10 u` after three sweeps (1342 of them above
/// `100 u`) and `6.7 u` after four, with none above `100 u`; five and more do not improve it and
/// cost orthogonality. `svd3_tests::four_sweeps_converge_and_three_do_not` pins both ends.
const SWEEPS: usize = 4;

/// `(c, s)` of the right rotation that orthogonalizes a column pair whose Gram entries are
/// `alpha = b_j . b_j`, `beta = b_k . b_k`, `gamma = b_j . b_k`.
///
/// Rutishauser's coefficients: `t` is the root of `t^2 + 2 zeta t - 1 = 0` of smaller magnitude,
/// `zeta = (alpha - beta) / (2 gamma)`, so `|theta| <= pi/4`. `sgn(+0) = +1`, so `alpha = beta`
/// gives `t = 1`, the 45-degree rotation the isotropic case needs.
///
/// Only `gamma = +-0` is guarded, where `zeta` would be `0/0` for `alpha = beta` and the pair is
/// already orthogonal. No negligibility threshold is typed (`0024` decision 2): a tiny `gamma`
/// yields a correct small `t`, or that same correct 45-degree rotation, so a threshold-skip would
/// only change results. `1 + zeta^2 >= 1` and the denominator is never zero, so no arm needs
/// `branch`.
#[inline]
fn givens<S: Real>(alpha: S, beta: S, gamma: S) -> (S, S) {
    let (zero, one) = (S::zero(), S::one());
    let flat = gamma.abs().le(zero);
    let zeta = (alpha - beta) / S::select(flat, one, S::lit(2.0) * gamma);
    let t = one.copysign(zeta) / (zeta.abs() + (one + zeta * zeta).sqrt());
    let t = S::select(flat, zero, t);
    let c = one / (one + t * t).sqrt();
    (c, t * c)
}

/// The pair `(x, y)` under `G = [[c, -s], [s, c]]`: `x' = c x + s y`, `y' = c y - s x`.
///
/// This pairing of the update with [`givens`]' `t` is normative (`0024` decision 2). The other
/// one rotates by `-theta` and leaves the off-diagonal at `2 cos(2 theta) gamma`: measured at
/// `1.0` after eight sweeps, against `4.8 u` here.
#[inline]
fn rotate<S: Real>(x: Vec3<S>, y: Vec3<S>, c: S, s: S) -> (Vec3<S>, Vec3<S>) {
    (x.scale(c) + y.scale(s), y.scale(c) - x.scale(s))
}

/// One rotation of the column pair, applied to `B` and to `V` alike.
#[inline]
fn orthogonalize<S: Real>(
    bx: Vec3<S>,
    by: Vec3<S>,
    vx: Vec3<S>,
    vy: Vec3<S>,
) -> (Vec3<S>, Vec3<S>, Vec3<S>, Vec3<S>) {
    let (c, s) = givens(bx.norm_sq(), by.norm_sq(), bx.dot(by));
    let (bx, by) = rotate(bx, by, c, s);
    let (vx, vy) = rotate(vx, vy, c, s);
    (bx, by, vx, vy)
}

/// A column of `B` with its norm and its column of `V`, the unit the sort moves.
type Col<S> = (S, Vec3<S>, Vec3<S>);

/// Descending compare-exchange, the permutation's parity carried as a `+-1` factor.
///
/// `A = sum_i b_i v_i^T` is invariant under a permutation applied to both, so the exchange itself
/// needs no compensation; the parity is spent once, at the end.
#[inline]
fn exchange<S: Real>(a: Col<S>, b: Col<S>, parity: S) -> (Col<S>, Col<S>, S) {
    let swap = a.0.lt(b.0);
    (
        Col::<S>::blend(swap, b, a),
        Col::<S>::blend(swap, a, b),
        parity * S::select(swap, -S::one(), S::one()),
    )
}

/// The singular value decomposition `A = U diag(sigma) V^T`.
///
/// Returns `(u, sigma, v)` with `u` and `v` orthonormal and `det u = det v = +1`,
/// `sigma.0[0] >= sigma.0[1] >= |sigma.0[2]|`, and **only `sigma.0[2]` ever negative**: it carries
/// the reflection of an input with `det a < 0`. So the nearest rotation to `a` in the Frobenius
/// norm is `u * v.transpose()` directly, with no sign to repair (`docs/NUMERICS.md` §3.4 projects
/// onto SO(3) this way).
///
/// Four sweeps of a cyclic one-sided Hestenes Jacobi over the column pairs `(0,1), (0,2), (1,2)`,
/// no convergence loop: each sweep applies Rutishauser's right rotation to the columns of
/// `B` (initially `a`) and of `V` (initially the identity), which drives `b_j . b_k` to zero.
/// Then `sigma_k = |b_k|`, sorted descending, and `u_0 = unit(b_0)`,
/// `u_1 = unit(b_1 - u_0 (u_0 . b_1))`, `u_2 = u_0 x u_1`. `A^T A` is never formed: forming it
/// costs the right singular vectors of a small close pair a factor `sigma_1 / sigma_2`, which is
/// why `0024` amended `docs/PHASE2.md` §6 away from McAdams.
///
/// `det u = +1` holds because `u_2` is a cross product, `det v = +1` because `V` is a product of
/// plane rotations, and an odd sort permutation is absorbed by negating `sigma_2` together with
/// the third column of `v`. No determinant is evaluated and no rank is decided by a threshold.
///
/// # Domain
///
/// Every input is legal and nothing is asserted (`docs/API.md` R6, as
/// [`chol`](crate::chol), [`solve_cubic`](crate::solve_cubic) and [`eig3`](crate::eig3)). **`u`
/// and `v` are finite, and their columns are ordered and normalized, for every input** — a
/// rank-deficient, non-finite or out-of-range one included, since every column is a normalization
/// or a cross product of two of them. The singular values are not finite where the input is not.
///
/// The Gram entries are squares that are not scaled, so **each column's squared norm must be a
/// normal number**: about `1e-150 < m < 1e150` (`f64`) and `1e-19 < m < 1e19` (`f32`) for entries
/// of magnitude `m`, as [`eig3`](crate::eig3). Inside that range the recorded bounds
/// (`svd3_tests`, a fit to the measurements, not a bar: `docs/PHASE2.md` §6 names no constant) are
///
/// - `|sigma_k - sigma_k*| <= 24 u |A|`, `|A| = sigma_1`; measured worst 6.0 `u |A|` (`f64`) and
///   7.6 (`f32`) over the mpmath fixture, 16 on a planted isotropic spectrum, which is the twelve
///   rotations' accumulated rounding;
/// - an angle of `u_k`, `v_k` of at most `24 u |A| / gap_k`, `gap_k` the distance from `sigma_k` to
///   the nearest other singular value in magnitude; measured 2.5 and 1.7;
/// - `|U^T U - I|`, `|V^T V - I|`, `|det - 1| <= 32 u` and a reconstruction
///   `|U diag(sigma) V^T - A| <= 32 u |A|`; measured 6.0 `u` and 5.0 `u |A|`.
///
/// **Outside that range only finiteness and the ordering survive.** A column norm that overflows
/// or goes subnormal leaves `zeta` saturated, so the rotation that would have separated two nearly
/// parallel columns is skipped and the projection building the second column cancels: the
/// orthogonality of `u` reaches about `sqrt(u)`. Scaling the input by its largest entry first
/// would remove the cliff for one rounding per entry, which is `0023` (draft) question 4 for
/// [`eig3`](crate::eig3) and is not answered here (`0024` decision 5).
///
/// **The singular values are accurate in `u |A|`, not relatively** — except on the column-graded
/// family `A = B diag(d)` with `B` well conditioned, where they come out to 1.8 `u sigma_k`
/// *relatively*, `sigma_3 / sigma_1 = 1e-32` included. That family is where the one-sided form
/// earns its keep (`0024`).
///
/// A `Dual` result differentiates the arm taken. `sigma_k` is a norm, so its derivative is NaN
/// where it is zero (`Dual::sqrt` at 0, decided in `0020`) — that is, for every rank-deficient
/// input; and a singular value is not differentiable there, so the NaN is the report.
///
/// # Example
///
/// ```
/// use helicoid_linalg::{svd3, Mat3, Vector};
///
/// // A reflection: det = -1, so the sign lands in sigma_3 and nowhere else.
/// let a = Mat3::from_rows([
///     Vector([2.0_f64, 0.0, 0.0]),
///     Vector([0.0, 0.0, 3.0]),
///     Vector([0.0, 1.0, 0.0]),
/// ]);
/// let (u, sigma, v) = svd3(&a);
/// assert!((sigma.0[0] - 3.0).abs() < 1e-15);
/// assert!((sigma.0[1] - 2.0).abs() < 1e-15);
/// assert!((sigma.0[2] + 1.0).abs() < 1e-15); // negative: det a < 0
///
/// // det u = det v = +1, so u v^T is the nearest rotation without a repair.
/// let det = |m: &Mat3<f64>| m.col(0).dot(m.col(1).cross(m.col(2)));
/// assert!((det(&u) - 1.0).abs() < 1e-14 && (det(&v) - 1.0).abs() < 1e-14);
/// ```
#[inline]
pub fn svd3<S: Real>(a: &Mat3<S>) -> (Mat3<S>, Vec3<S>, Mat3<S>) {
    with_sweeps::<S, SWEEPS>(a)
}

/// [`svd3`] with the sweep count open, so `svd3_tests` can measure that four sweeps converge and
/// three do not (`0024` decision 3). Not public: `N` is not a parameter of the decomposition.
pub(crate) fn with_sweeps<S: Real, const N: usize>(a: &Mat3<S>) -> (Mat3<S>, Vec3<S>, Mat3<S>) {
    let (zero, one) = (S::zero(), S::one());
    let (mut b0, mut b1, mut b2) = (a.col(0), a.col(1), a.col(2));
    let mut v0 = Vector([one, zero, zero]);
    let mut v1 = Vector([zero, one, zero]);
    let mut v2 = Vector([zero, zero, one]);

    for _ in 0..N {
        (b0, b1, v0, v1) = orthogonalize(b0, b1, v0, v1);
        (b0, b2, v0, v2) = orthogonalize(b0, b2, v0, v2);
        (b1, b2, v1, v2) = orthogonalize(b1, b2, v1, v2);
    }

    // Descending by column norm, so `|sigma_3|` is the smallest: the ordering `u v^T` needs to be
    // the nearest rotation.
    let (c0, c1, c2) = (
        (b0.norm(), b0, v0),
        (b1.norm(), b1, v1),
        (b2.norm(), b2, v2),
    );
    let (c0, c1, parity) = exchange(c0, c1, one);
    let (c1, c2, parity) = exchange(c1, c2, parity);
    let (c0, c1, parity) = exchange(c0, c1, parity);

    // Both frames are built the same way, from the two leading columns and a cross product: the
    // third column is then `+-` the one it replaces, so `det = +1` holds by construction and no
    // input can break it. The fallbacks of `unit_or` only have to be orthonormal, which
    // `any_orthogonal` is by construction, so `U` and `V` are finite even where the sweeps
    // produced nothing usable — a non-finite entry poisons `V`'s accumulated rotations.
    let frame = |x: Vec3<S>, y: Vec3<S>| {
        let a0 = unit_or(x, Vector([one, zero, zero]));
        let fallback = any_orthogonal(a0);
        let a1 = unit_or(y - a0.scale(a0.dot(y)), fallback);
        [a0, a1, a0.cross(a1)]
    };
    let [u0, u1, u2] = frame(c0.1, c1.1);
    let [v0, v1, v2] = frame(c0.2, c1.2);

    // `V` was a product of plane rotations, so `v_0 x v_1` is its own third column times the
    // sort's parity; `b_2 = sigma_3 u_2` up to the orthogonality residual, so the sign of
    // `u_2 . b_2` is the reflection of the input. Their product leaves `sigma_3 v_3^T` unchanged
    // and puts the whole sign in `sigma_3`, where §6 wants it.
    let sigma2 = c2.0.copysign(u2.dot(c2.1) * parity);
    (
        Mat3::from_cols([u0, u1, u2]),
        Vector([c0.0, c1.0, sigma2]),
        Mat3::from_cols([v0, v1, v2]),
    )
}
