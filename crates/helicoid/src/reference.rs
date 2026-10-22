//! The reference twins: the obvious, slow implementations that the fast composites are
//! proptested against (`NUMERICS.md` §14, D6).
//!
//! A twin is the definition of *correct* for its fast counterpart, so it is public and it is not
//! optimized. Each writes its dense result into caller memory (`docs/API.md` R5: no function
//! returns a dense matrix) and keeps its scratch on the stack, sized by a const parameter `D` that
//! must equal the dimension, asserted at monomorphization (`cargo build` and `cargo test`, not
//! `cargo check`), as for [`Jac::sandwich`]:
//!
//! ```compile_fail,E0080
//! use helicoid::{reference, Jac, SEn3Jac};
//! use helicoid_linalg::StridedMut;
//! let j: SEn3Jac<f64, 1> = Jac::identity();
//! let mut out = [0.0; 49];
//! reference::sen3jac_mul::<f64, 1, 7>(&j, &j, &mut StridedMut::col_major(&mut out, 7, 7));
//! ```
//!
//! Positive control:
//!
//! ```
//! use helicoid::{reference, Jac, SEn3Jac};
//! use helicoid_linalg::StridedMut;
//! let j: SEn3Jac<f64, 1> = Jac::identity();
//! let mut out = [0.0; 36];
//! reference::sen3jac_mul::<f64, 1, 6>(&j, &j, &mut StridedMut::col_major(&mut out, 6, 6));
//! ```

use crate::dualmat::{dof, SEn3Jac};
use crate::product::ProductJac;
use crate::side::Side;
use crate::so3::SO3Tangent;
use crate::traits::{Jac, LieGroup, Tangent};
use crate::{Chart, Quat, RnTangent, S2Chart, Sim3Jac, S2, SE3, SO3};
use helicoid_linalg::{Mask, Mat3, Matrix, Point3, Precision, Real, StridedMut, Vec3, Vector};

/// The terms added left to right from the first, `+0` when there are none.
///
/// The same association as `helicoid_linalg`'s `vector::sum`, which is `pub(crate)` there and so out
/// of reach from here. `Matrix`'s `Mul` sums through that one and every twin below sums through this
/// one, and `integer_products_and_inverses_are_exact`, `sandwich_of_integers_is_exact` and
/// `apply_and_apply_transpose_are_the_dense_matrix_and_its_transpose` all assert the two agree bit
/// for bit; a signed zero tells the seeds apart. So the two must change together, and neither may
/// change alone. Making linalg's public instead is a new public item and owes a record
/// (`docs/API.md` §6).
fn sum<S: Real>(mut terms: impl Iterator<Item = S>) -> S {
    match terms.next() {
        Some(first) => terms.fold(first, |acc, t| acc + t),
        None => S::zero(),
    }
}

/// `Γ_m(φ) v = Σ Wⁿ v/(n+m)!`, the definition (`NUMERICS.md` §7), summed left to right over its
/// first [`GAMMA_TERMS`] terms with `Wⁿ v` as repeated cross products: the twin of
/// [`so3::gamma_apply_jacobians`](crate::so3::gamma_apply_jacobians), run on `Dual<S, 3>` for
/// its Jacobian (§14).
///
/// # Domain
///
/// `m!` exact at `S`, as `S::lit` requires: `m ≤ 18` at binary64 and `m ≤ 13` at binary32,
/// `debug_assert!`ed (`0066`). `θ ≤ π` for the stated accuracy: the last term is below
/// `π⁴⁰/40! < 10⁻²⁷` of `v` there. Beyond it the truncation is not bounded by this constant, and
/// the terms grow to `θⁿ/n!` before they fall, which cancels.
pub fn gamma_apply<S: Real>(m: usize, phi: &SO3Tangent<S>, v: Vec3<S>) -> Vec3<S> {
    debug_assert!(
        m <= match S::PRECISION {
            Precision::F64 => 18,
            Precision::F32 => 13,
        },
        "reference::gamma_apply: m! is not exact at this precision"
    );
    let div = |w: Vec3<S>, k: f64| Vector(w.0.map(|x| x / S::lit(k)));
    let mut term = div(v, factorial(m));
    let mut acc = term;
    for n in 1..GAMMA_TERMS {
        term = div(phi.phi.cross(term), (n + m) as f64);
        acc = acc + term;
    }
    acc
}

/// The number of terms [`gamma_apply`] sums.
pub const GAMMA_TERMS: usize = 40;

/// `m!` for the `m` of [`gamma_apply`], exact in binary64 through `m = 18` and binary32 through 13.
fn factorial(m: usize) -> f64 {
    (1..=m).fold(1.0, |f, k| f * k as f64)
}

/// The group geodesic, the definition: `x0 ⊕_R ((x1 ⊖_R x0) · t)` (`NUMERICS.md` §10).
///
/// `LieGroup::geodesic`'s provided body is a call to this, so the two are one expression
/// (`PHASE4.md` §1.1) and a group's fast twin has something to be compared against that its own
/// override cannot shadow. The only twin here that is not dense: a geodesic is a group element,
/// so `API.md` R5 does not apply.
///
/// # Domain
///
/// [`LieGroup::geodesic`]'s.
pub fn geodesic<S: Real, G: LieGroup<S>>(x0: &G, x1: &G, t: S) -> G {
    x0.rplus(&x1.rminus(x0).scale(t))
}

/// `SO3::inv_act_jacobians`' twin (`NUMERICS.md` §14, `0103`): the composition it replaces,
/// `act_jacobians(R⁻¹, v) · inverse_jacobian(R)`, each entry summed left to right from its first
/// term. Writes `∂(R⁻¹ v)/∂R` in the convention of `Sd` to the `3 x 3` view `out`. `∂/∂v` has no
/// twin: `act_jacobians(R⁻¹, v)`'s `R(q*)` is `R(q)ᵀ` to the bit.
///
/// # Domain
///
/// `out` is `3 x 3`, checked by `debug_assert!`; a smaller view panics in the strided access (D11).
pub fn so3_inv_act_jacobian<S: Real, Sd: Side>(
    r: &SO3<S>,
    v: Vec3<S>,
    out: &mut StridedMut<'_, S>,
) {
    debug_assert!(
        out.rows() == 3 && out.cols() == 3,
        "so3_inv_act_jacobian: the view is not 3 x 3"
    );
    let a = r.inverse().act_jacobians::<Sd>(v).0;
    let b = dense::<S, _, Mat3<S>, 3>(&r.inverse_jacobian::<Sd>());
    product_into(&a, &b, out);
}

/// `SE3::inv_act_jacobians`' twin (`NUMERICS.md` §14, `0103`): `act_jacobians(X⁻¹, p) ·
/// inverse_jacobian(X)`, the second written dense, each entry of the `3 x 6` product summed left
/// to right from its first term. Writes `∂(X⁻¹ p)/∂X` in the convention of `Sd` to the `3 x 6`
/// view `out`; `∂/∂p` has no twin, for [`so3_inv_act_jacobian`]'s reason.
///
/// # Domain
///
/// `out` is `3 x 6`, checked by `debug_assert!`; a smaller view panics in the strided access (D11).
pub fn se3_inv_act_jacobian<S: Real, Sd: Side>(
    x: &SE3<S>,
    p: Point3<S>,
    out: &mut StridedMut<'_, S>,
) {
    debug_assert!(
        out.rows() == 3 && out.cols() == 6,
        "se3_inv_act_jacobian: the view is not 3 x 6"
    );
    let a = x.inverse().act_jacobians::<Sd>(p).0;
    let b = dense::<S, _, SEn3Jac<S, 1>, 6>(&x.inverse_jacobian::<Sd>());
    product_into(&a, &b, out);
}

/// `a · b` for a `3 x D` and a row-major `D x D`, into `out`, summed left to right.
fn product_into<S: Real, const D: usize>(
    a: &Matrix<S, 3, D>,
    b: &[[S; D]; D],
    out: &mut StridedMut<'_, S>,
) {
    for r in 0..3 {
        for c in 0..D {
            out.set(
                r,
                c,
                sum(b.iter().enumerate().map(|(k, bk)| a.get(r, k) * bk[c])),
            );
        }
    }
}

/// `S2Chart::local`'s twin, its definition (`NUMERICS.md` §14, `0063` decision 4): the minimal
/// rotation taking `n` to `m` as the quaternion `(1 + n·m, n × m)`, normalized, its `Log` through
/// `SO3::log`'s `atan2`, projected on the chart's basis. Not the fast path's `α/s` from
/// `log_ratio`, nor its `Bᵀ(α/s · n × m)`: this one halves the angle, which is `r(s², 1 + w)`, the
/// form `0063` decision 3 rejects as 3 to 5 times less accurate near the antipode (CH.9(c)).
///
/// # Domain
///
/// `S2Chart::local`'s, and `(1 + n·m)² + ‖n × m‖²` normal, for the normalization.
pub fn s2_local<S: Real>(chart: &S2Chart<S>, m: &S2<S>) -> RnTangent<S, 2> {
    let (n, m) = (chart.base().vec(), m.vec());
    let [x, y, z] = n.cross(m).0;
    let q = Quat::from_wxyz_normalized(S::one() + n.dot(m), x, y, z);
    let phi = SO3::from_quat_unchecked(q).log().phi;
    RnTangent {
        rho: Vector([chart.b1.dot(phi), chart.b2.dot(phi)]),
    }
}

/// The dense image of `j` in a `D x D` scratch, through the public [`Jac::write_dense`].
///
/// The scratch is NaN-poisoned, so an entry `write_dense` leaves unwritten reaches the result
/// instead of passing as a structural zero. `(rs, cs)` is the orientation, and the only thing the
/// two callers differ in: [`dense`] passes the `(D, 1)` of `StridedMut::row_major`, so `m[r][c]` is
/// entry `(r, c)`, and the hand cases of `laws` pass the `(1, D)` of `col_major` for `m[c][r]`. The
/// strides and not a constructor, because `StridedMut`'s lifetime is a parameter of the type, so
/// `row_major` is early-bound and cannot be passed as a `for<'a>` callback. One body, so the poison
/// has one place to be kept rather than two to be kept in step.
pub(crate) fn dense_oriented<S: Real, T: Tangent<S>, J: Jac<S, T>, const D: usize>(
    j: &J,
    rs: usize,
    cs: usize,
) -> [[S; D]; D] {
    let mut m = [[S::zero() / S::zero(); D]; D];
    j.write_dense(&mut StridedMut::with_strides(
        m.as_flattened_mut(),
        D,
        D,
        rs,
        cs,
    ));
    m
}

/// [`dense_oriented`] row-major (`m[r][c]`), the orientation every twin here reads.
pub(crate) fn dense<S: Real, T: Tangent<S>, J: Jac<S, T>, const D: usize>(j: &J) -> [[S; D]; D] {
    dense_oriented(j, D, 1)
}

/// The dense `M Σ Mᵀ` of the row-major `m`, written to the `D x D` view `out`.
///
/// Entry `(i, k)` is `Σ_p Σ_q m_ip Σ_pq m_kq`: the inner sum over `q` left to right from its first
/// term, then those `D` sums over `p` the same way. Every one of the `D⁴` products is formed,
/// structural zeros included; `Σ` need not be symmetric and the result is not symmetrized. Every
/// `sandwich` twin is this applied to its own dense image.
fn sandwich_dense<S: Real, const D: usize>(
    m: &[[S; D]; D],
    cov: &Matrix<S, D, D>,
    out: &mut StridedMut<'_, S>,
) {
    for (i, mi) in m.iter().enumerate() {
        for (k, mk) in m.iter().enumerate() {
            let row = |p: usize| sum((0..D).map(|q| mi[p] * cov.get(p, q) * mk[q]));
            out.set(i, k, sum((0..D).map(row)));
        }
    }
}

/// The twin of [`SEn3Jac::mul`](Jac::mul): the dense product `dense(a) * dense(b)`, every one of
/// the `D³` products formed, structural zeros included, each entry summed left to right from
/// its first term. Writes it to the `D x D` view `out`.
///
/// # Domain
///
/// `D == 3 + 3N` (a build-time assertion) and `out` is `D x D`, checked by `debug_assert!`; a
/// smaller view panics in the strided access, the one documented panic class (D11).
pub fn sen3jac_mul<S: Real, const N: usize, const D: usize>(
    a: &SEn3Jac<S, N>,
    b: &SEn3Jac<S, N>,
    out: &mut StridedMut<'_, S>,
) {
    const { assert!(D == dof::<S, N>()) };
    debug_assert!(
        out.rows() == D && out.cols() == D,
        "sen3jac_mul: the view is not D x D"
    );
    mul_dense::<S, _, _, D>(a, b, out);
}

/// The dense product of two `Jac`s, row by column, summed left to right: the body of every
/// `*jac_mul` twin.
fn mul_dense<S: Real, T: Tangent<S>, J: Jac<S, T>, const D: usize>(
    a: &J,
    b: &J,
    out: &mut StridedMut<'_, S>,
) {
    let (x, y) = (dense::<S, _, _, D>(a), dense::<S, _, _, D>(b));
    for (r, row) in x.iter().enumerate() {
        for c in 0..D {
            out.set(r, c, sum(row.iter().zip(&y).map(|(a, yk)| *a * yk[c])));
        }
    }
}

/// The twin of [`SEn3Jac::inverse`](Jac::inverse): Gauss–Jordan elimination of `dense(j)`
/// beside the identity, which becomes the inverse. Writes it to the `D x D` view `out`.
///
/// The pivot of column `k` is the entry of largest magnitude in rows `k..`, found by comparing
/// row `k` with each later row in turn and exchanging the two where the later pivot is larger
/// (`S::select`, so no comparison reaches an `if`, `docs/API.md` R4). Its forward error is
/// `O(κ u)` for the conditioning `κ = ‖M‖_F ‖M⁻¹‖_F` of the *dense* matrix `M`, so it and the
/// structured inverse agree to a few `κ u`, not to a few `u`: the value differed by at most
/// `1.25 κ ‖M⁻¹‖_F u` in `‖·‖_F` over `10^6` seeded cases per scalar and `N = 1, 2, 3`, with
/// entries uniform on `[-1, 1)` and `κ u <= 1e-3` (`sen3jac_inverse_matches_reference`, which
/// rejects the other inputs). A nearly singular `A` is not covered by that measurement.
///
/// # Domain
///
/// `D == 3 + 3N` (a build-time assertion), `out` is `D x D` and every computed pivot is nonzero;
/// the last two are checked by `debug_assert!` (a NaN fails the pivot check). A release build
/// divides by a zero pivot and returns non-finite entries. A smaller `out` panics in the strided
/// access (D11).
///
/// ```compile_fail,E0080
/// use helicoid::{reference, Jac, SEn3Jac};
/// use helicoid_linalg::StridedMut;
/// let j: SEn3Jac<f64, 1> = Jac::identity();
/// let mut out = [0.0; 49];
/// reference::sen3jac_inverse::<f64, 1, 7>(&j, &mut StridedMut::col_major(&mut out, 7, 7));
/// ```
///
/// Positive control:
///
/// ```
/// use helicoid::{reference, Jac, SEn3Jac};
/// use helicoid_linalg::StridedMut;
/// let j: SEn3Jac<f64, 1> = Jac::identity();
/// let mut out = [0.0; 36];
/// reference::sen3jac_inverse::<f64, 1, 6>(&j, &mut StridedMut::col_major(&mut out, 6, 6));
/// ```
pub fn sen3jac_inverse<S: Real, const N: usize, const D: usize>(
    j: &SEn3Jac<S, N>,
    out: &mut StridedMut<'_, S>,
) {
    const { assert!(D == dof::<S, N>()) };
    debug_assert!(
        out.rows() == D && out.cols() == D,
        "sen3jac_inverse: the view is not D x D"
    );
    inverse_dense::<S, _, _, D>(j, out);
}

/// Gauss–Jordan with partial pivoting on the dense matrix: the body of every `*jac_inverse` twin.
fn inverse_dense<S: Real, T: Tangent<S>, J: Jac<S, T>, const D: usize>(
    j: &J,
    out: &mut StridedMut<'_, S>,
) {
    let mut m = dense::<S, _, _, D>(j);
    for r in 0..D {
        for c in 0..D {
            out.set(r, c, if r == c { S::one() } else { S::zero() });
        }
    }
    for k in 0..D {
        for r in k + 1..D {
            let swap = m[k][k].abs().lt(m[r][k].abs());
            let (head, tail) = m.split_at_mut(r);
            for (p, q) in head[k].iter_mut().zip(tail[0].iter_mut()) {
                (*p, *q) = (S::select(swap, *q, *p), S::select(swap, *p, *q));
            }
            for c in 0..D {
                let (p, q) = (out.get(k, c), out.get(r, c));
                out.set(k, c, S::select(swap, q, p));
                out.set(r, c, S::select(swap, p, q));
            }
        }
        let pivot = m[k][k];
        debug_assert!(
            S::zero().lt(pivot.abs()).all(),
            "reference: a zero or NaN pivot"
        );
        for x in &mut m[k] {
            *x = *x / pivot;
        }
        for c in 0..D {
            out.set(k, c, out.get(k, c) / pivot);
        }
        let row = m[k];
        for (r, mr) in m.iter_mut().enumerate().filter(|&(r, _)| r != k) {
            let f = mr[k];
            for (x, y) in mr.iter_mut().zip(&row) {
                *x = *x - f * *y;
            }
            for c in 0..D {
                out.set(r, c, out.get(r, c) - f * out.get(k, c));
            }
        }
    }
}

/// The twin of [`SEn3Jac::apply`](Jac::apply): the dense product `dense(j) x`, every one of the
/// `D²` products formed, structural zeros included, each entry summed left to right from its first
/// term. Writes it to `out`.
///
/// The twin takes and returns the dense components, not a [`SEn3Tangent`](crate::SEn3Tangent): the
/// definition of *correct* here is the matrix-vector product, and the structure of the operand is
/// what the fast path is allowed to exploit. A caller compares through
/// [`write_dense`](crate::Tangent::write_dense).
///
/// # Domain
///
/// `D == 3 + 3N` (a build-time assertion) and `x.len() == out.len() == D`, checked by
/// `debug_assert!`. A release build never panics: a short `x` drops the products past its end and
/// a wrongly sized `out` is filled to `min(out.len(), D)`.
pub fn sen3jac_apply<S: Real, const N: usize, const D: usize>(
    j: &SEn3Jac<S, N>,
    x: &[S],
    out: &mut [S],
) {
    const { assert!(D == dof::<S, N>()) };
    debug_assert!(
        x.len() == D && out.len() == D,
        "sen3jac_apply: an operand is not D long"
    );
    apply_dense::<S, _, _, D>(j, x, out);
}

/// `J x` on the dense matrix, each row summed left to right.
fn apply_dense<S: Real, T: Tangent<S>, J: Jac<S, T>, const D: usize>(
    j: &J,
    x: &[S],
    out: &mut [S],
) {
    let m = dense::<S, _, _, D>(j);
    for (o, row) in out.iter_mut().zip(&m) {
        *o = sum(row.iter().zip(x).map(|(a, b)| *a * *b));
    }
}

/// The twin of [`SEn3Jac::apply_transpose`](Jac::apply_transpose): the dense product
/// `dense(j)ᵀ x`, formed and summed as [`sen3jac_apply`] does, down the columns instead of along
/// the rows. Writes it to `out`.
///
/// The fast path does not associate this sum the same way — it adds `Aᵀ φ` and then one whole
/// `B_iᵀ ρ_i` per block, where a dense row is one flat sum of `D` terms — so the two agree to a
/// few `u`, not bit for bit. `Tangent::dot_acc`'s index order is normative (`0025`); the
/// association inside a block product is not.
///
/// # Domain
///
/// As [`sen3jac_apply`].
pub fn sen3jac_apply_transpose<S: Real, const N: usize, const D: usize>(
    j: &SEn3Jac<S, N>,
    x: &[S],
    out: &mut [S],
) {
    const { assert!(D == dof::<S, N>()) };
    debug_assert!(
        x.len() == D && out.len() == D,
        "sen3jac_apply_transpose: an operand is not D long"
    );
    apply_transpose_dense::<S, _, _, D>(j, x, out);
}

/// `Jᵀ x` on the dense matrix, each column summed top to bottom.
fn apply_transpose_dense<S: Real, T: Tangent<S>, J: Jac<S, T>, const D: usize>(
    j: &J,
    x: &[S],
    out: &mut [S],
) {
    let m = dense::<S, _, _, D>(j);
    for (r, o) in (0..D).zip(out.iter_mut()) {
        *o = sum(m.iter().zip(x).map(|(row, b)| row[r] * *b));
    }
}

/// The twin of [`SEn3Jac::sandwich`](Jac::sandwich): the dense `dense(j) Σ dense(j)ᵀ`, every one
/// of the `D⁴` products formed, structural zeros included. Writes it to the `D x D` view `out`.
///
/// Entry `(i, k)` is `Σ_p Σ_q J_ip Σ_pq J_kq`: the inner sum over `q` left to right from its first
/// term, then those `D` sums over `p` the same way. `Σ` need not be symmetric and the result is
/// not symmetrized, as for the fast path. This is also the twin `NUMERICS.md` §14 owes
/// `Gaussian::to_left`/`to_right` at Phase 5 — `Ad Σ Adᵀ` is this product with `j = Ad`.
///
/// # Domain
///
/// `D == 3 + 3N` (a build-time assertion) and `out` is `D x D`, checked by `debug_assert!`; a
/// smaller view panics in the strided access, the one documented panic class (D11). Here `D` ties
/// three things — the scratch, the type of `cov` and the view — so the assertion is pinned:
///
/// ```compile_fail,E0080
/// use helicoid::{reference, Jac, SEn3Jac};
/// use helicoid_linalg::{Matrix, StridedMut};
/// let j: SEn3Jac<f64, 1> = Jac::identity();
/// let mut out = [0.0; 49];
/// let cov = Matrix::<f64, 7, 7>::identity();
/// reference::sen3jac_sandwich::<f64, 1, 7>(&j, &cov, &mut StridedMut::col_major(&mut out, 7, 7));
/// ```
///
/// Positive control:
///
/// ```
/// use helicoid::{reference, Jac, SEn3Jac};
/// use helicoid_linalg::{Matrix, StridedMut};
/// let j: SEn3Jac<f64, 1> = Jac::identity();
/// let mut out = [0.0; 36];
/// let cov = Matrix::<f64, 6, 6>::identity();
/// reference::sen3jac_sandwich::<f64, 1, 6>(&j, &cov, &mut StridedMut::col_major(&mut out, 6, 6));
/// assert_eq!(out[0].to_bits(), 1.0_f64.to_bits());
/// ```
pub fn sen3jac_sandwich<S: Real, const N: usize, const D: usize>(
    j: &SEn3Jac<S, N>,
    cov: &Matrix<S, D, D>,
    out: &mut StridedMut<'_, S>,
) {
    const { assert!(D == dof::<S, N>()) };
    debug_assert!(
        out.rows() == D && out.cols() == D,
        "sen3jac_sandwich: the view is not D x D"
    );
    sandwich_dense(&dense::<S, _, _, D>(j), cov, out);
}

/// The twin of [`ProductJac::sandwich`](Jac::sandwich): the dense `dense(j) Σ dense(j)ᵀ`, every
/// one of the `D⁴` products formed and summed as [`sen3jac_sandwich`] forms and sums them. Writes
/// it to the `D x D` view `out`.
///
/// The fast path never builds a dense `J`: it applies the factors' blocks to the columns of `Σ`
/// and then to the rows of `J Σ`. It therefore associates its sums differently — a block product
/// against one flat row of `D` terms — so the two agree to a few `u`, not to the bit.
///
/// # Domain
///
/// `D == DOF` (a build-time assertion) and `out` is `D x D`, checked by `debug_assert!`; a smaller
/// view panics in the strided access, the one documented panic class (D11):
///
/// ```compile_fail,E0080
/// use helicoid::{reference, Jac, ProductJac, RnJac};
/// use helicoid_linalg::{Matrix, StridedMut};
/// let j: ProductJac<RnJac<f64, 2>, RnJac<f64, 3>> = Jac::identity();
/// let (mut out, cov) = ([0.0; 36], Matrix::<f64, 6, 6>::identity());
/// let mut view = StridedMut::col_major(&mut out, 6, 6);
/// reference::productjac_sandwich::<f64, _, _, _, _, 6>(&j, &cov, &mut view);
/// ```
///
/// Positive control:
///
/// ```
/// use helicoid::{reference, Jac, ProductJac, RnJac};
/// use helicoid_linalg::{Matrix, StridedMut};
/// let j: ProductJac<RnJac<f64, 2>, RnJac<f64, 3>> = Jac::identity();
/// let (mut out, cov) = ([0.0; 25], Matrix::<f64, 5, 5>::identity());
/// let mut view = StridedMut::col_major(&mut out, 5, 5);
/// reference::productjac_sandwich::<f64, _, _, _, _, 5>(&j, &cov, &mut view);
/// assert_eq!(out[0].to_bits(), 1.0_f64.to_bits());
/// ```
pub fn productjac_sandwich<S, TA, TB, JA, JB, const D: usize>(
    j: &ProductJac<JA, JB>,
    cov: &Matrix<S, D, D>,
    out: &mut StridedMut<'_, S>,
) where
    S: Real,
    TA: Tangent<S>,
    TB: Tangent<S>,
    JA: Jac<S, TA>,
    JB: Jac<S, TB>,
{
    const { assert!(D == <(TA, TB) as Tangent<S>>::DOF) };
    debug_assert!(
        out.rows() == D && out.cols() == D,
        "productjac_sandwich: the view is not D x D"
    );
    sandwich_dense(&dense::<S, (TA, TB), ProductJac<JA, JB>, D>(j), cov, out);
}

/// [`sen3jac_mul`]'s twin for [`Sim3Jac`]: the dense `7 x 7` product
/// (`NUMERICS.md` §14, `0070`).
///
/// # Domain
///
/// `out` is `7 x 7`, checked by `debug_assert!`.
pub fn sim3jac_mul<S: Real>(a: &Sim3Jac<S>, b: &Sim3Jac<S>, out: &mut StridedMut<'_, S>) {
    debug_assert!(
        out.rows() == 7 && out.cols() == 7,
        "sim3jac_mul: the view is not 7 x 7"
    );
    mul_dense::<S, _, _, 7>(a, b, out);
}

/// The dense Gauss–Jordan inverse of a [`Sim3Jac`], the twin of
/// `Sim3Jac::inverse` and of `Sim3::jr_inv`'s block algebra (`NUMERICS.md` §14, `0070`).
///
/// # Domain
///
/// As [`sen3jac_inverse`].
pub fn sim3jac_inverse<S: Real>(j: &Sim3Jac<S>, out: &mut StridedMut<'_, S>) {
    debug_assert!(
        out.rows() == 7 && out.cols() == 7,
        "sim3jac_inverse: the view is not 7 x 7"
    );
    inverse_dense::<S, _, _, 7>(j, out);
}

/// The dense `J x` of a [`Sim3Jac`] (`NUMERICS.md` §14).
///
/// # Domain
///
/// `x.len() == out.len() == 7`, checked by `debug_assert!`.
pub fn sim3jac_apply<S: Real>(j: &Sim3Jac<S>, x: &[S], out: &mut [S]) {
    debug_assert!(
        x.len() == 7 && out.len() == 7,
        "sim3jac_apply: an operand is not 7 long"
    );
    apply_dense::<S, _, _, 7>(j, x, out);
}

/// The dense `Jᵀ x` of a [`Sim3Jac`] (`NUMERICS.md` §14).
///
/// # Domain
///
/// As [`sim3jac_apply`].
pub fn sim3jac_apply_transpose<S: Real>(j: &Sim3Jac<S>, x: &[S], out: &mut [S]) {
    debug_assert!(
        x.len() == 7 && out.len() == 7,
        "sim3jac_apply_transpose: an operand is not 7 long"
    );
    apply_transpose_dense::<S, _, _, 7>(j, x, out);
}

/// The dense `J Σ Jᵀ` of a [`Sim3Jac`] (`NUMERICS.md` §14).
///
/// # Domain
///
/// `out` is `7 x 7`, checked by `debug_assert!`.
pub fn sim3jac_sandwich<S: Real>(
    j: &Sim3Jac<S>,
    cov: &Matrix<S, 7, 7>,
    out: &mut StridedMut<'_, S>,
) {
    debug_assert!(
        out.rows() == 7 && out.cols() == 7,
        "sim3jac_sandwich: the view is not 7 x 7"
    );
    sandwich_dense(&dense::<S, _, _, 7>(j), cov, out);
}
