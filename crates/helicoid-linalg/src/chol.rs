//! Fixed-size Cholesky factorization and triangular solves (`docs/PHASE2.md` §4).
//!
//! Every inner product is accumulated left to right from the first term ([`crate::vector`]'s
//! reduction), so a result is a fixed sequence of `S` operations (D16). Nothing here branches on
//! a float: the comparisons (the pivot test, the finiteness of an entry) are masks.

use crate::matrix::Matrix;
use crate::real::{is_finite, Mask, Real};
use crate::vector::{sum, Vector};

/// The Cholesky factor of `a` and whether `a` is positive definite.
///
/// Returns the lower-triangular `L` (exact `+0` above the diagonal) with `L L^T = a`, and a mask
/// that is set iff every pivot was strictly positive. Only the lower triangle of `a` is read; the
/// caller owns symmetry. The order is Cholesky-Crout, column by column:
///
/// ```text
/// d_j    = a_jj - sum_{k<j} L_jk^2          pivot, the sums run in k order
/// L_jj   = sqrt(d_j)
/// L_ij   = (a_ij - sum_{k<j} L_ik L_jk) / L_jj      i > j
/// ```
///
/// The pivot test is `0 < d_j` and `d_j` finite, so it is false for a zero pivot (`+0` or `-0`: a
/// singular matrix is not positive definite), for a negative one, for a NaN and for `+inf`, and
/// true for the smallest subnormal. An entry `L_ij` that is not finite (an overflowing quotient,
/// or a NaN or `inf` read from `a`) also clears the mask and is stored as `+0`. The mask is a
/// statement about the *computed* values, not a certificate about `a`. When it is set,
/// `L L^T = a + E` with `|E| <= gamma_{N+1} |L| |L|^T` entrywise (Higham, *Accuracy and
/// Stability of Numerical Algorithms*, Thm 10.3), however ill-conditioned `a` is: the factor is
/// backward stable. A matrix that is positive semidefinite in exact arithmetic (rank-deficient)
/// or positive definite with `cond(a)` near `1/u` can come back either way, because a pivot that
/// is zero or tiny in exact arithmetic is rounding noise; a clear mask means "a pivot was not
/// positive as computed, or a factor entry overflowed". The derivation and the bounds are `NUMERICS.md` §15.
///
/// Every lane computes every arm, so nothing is fed a value it cannot take. `sqrt` gets
/// `select(ok, d, 1)`, never a negative or NaN argument (a `Dual` derivative stays finite), and a
/// failed pivot gives `L_jj = 1` and zeros below it in its column, so garbage does not feed later
/// columns (without the zeros it squares per column and overflows `f32` at `N = 6` on entries
/// below 20).
///
/// # Domain
///
/// Positive definite is the domain on which `L` is the Cholesky factor (`NUMERICS.md` §12). Any
/// `a` is legal input and nothing is asserted: the mask is the report (`docs/API.md` R4, R6).
/// Outside the domain, in release, the result is `(L, mask = false)` with unspecified values, and
/// `L` is finite (its value lane, for a `Dual`) for every `a`, NaN and `inf` in the lower
/// triangle included: no non-finite entry is ever stored. Values above the diagonal are exact
/// `+0`. The finiteness of a `Dual`'s derivative lanes is that of `Real::sqrt` and the quotient
/// (`NUMERICS.md` §12).
///
/// # Example
///
/// ```
/// use helicoid_linalg::{chol, Mat2, Vector};
///
/// let a = Mat2::from_cols([Vector([4.0_f64, 2.0]), Vector([2.0, 5.0])]);
/// let (l, pd) = chol(&a);
/// assert!(pd && l.get(0, 0) == 2.0 && l.get(1, 0) == 1.0 && l.get(1, 1) == 2.0);
/// ```
#[inline]
pub fn chol<S: Real, const N: usize>(a: &Matrix<S, N, N>) -> (Matrix<S, N, N>, S::Mask) {
    let mut cols = [[S::zero(); N]; N];
    let mut pd = S::zero().le(S::zero());
    for j in 0..N {
        let d = a.get(j, j) - sum((0..j).map(|k| cols[k][j] * cols[k][j]));
        let ok = S::zero().lt(d).and(is_finite(d));
        pd = pd.and(ok);
        let ljj = S::select(ok, d, S::one()).sqrt();
        cols[j][j] = ljj;
        for i in j + 1..N {
            let s = sum((0..j).map(|k| cols[k][i] * cols[k][j]));
            let x = (a.get(i, j) - s) / ljj;
            let fin = is_finite(x);
            pd = pd.and(fin.or(ok.not()));
            cols[j][i] = S::select(ok.and(fin), x, S::zero());
        }
    }
    (Matrix::from_cols(cols.map(Vector)), pd)
}

/// Solves `l x = b` by forward substitution, reading only the lower triangle of `l`.
///
/// `x_i = (b_i - sum_{k<i} l_ik x_k) / l_ii`, the sum in `k` order. Backward stable:
/// `(l + E) x = b` with `|E| <= gamma_N |l|` (Higham, Thm 8.5).
///
/// # Domain
///
/// Every `l_ii` is nonzero and not NaN, checked by `debug_assert!`. A release build never panics:
/// a zero diagonal divides by zero and the result is `inf` or NaN.
#[inline]
pub fn solve_lower<S: Real, const N: usize>(l: &Matrix<S, N, N>, b: Vector<S, N>) -> Vector<S, N> {
    let mut x = [S::zero(); N];
    for i in 0..N {
        let lii = l.get(i, i);
        debug_assert!(lii.abs().value_f64() > 0.0, "solve_lower: zero diagonal");
        x[i] = (b.0[i] - sum((0..i).map(|k| l.get(i, k) * x[k]))) / lii;
    }
    Vector(x)
}

/// Solves `u x = b` by back substitution, reading only the upper triangle of `u`. For `A x = b`
/// from a Cholesky factor use [`chol_solve`], which does not build `l.transpose()`.
///
/// `x_i = (b_i - sum_{k>i} u_ik x_k) / u_ii`, the sum in increasing `k`, `i` from `N - 1` down.
/// Backward stable like [`solve_lower`].
///
/// # Domain
///
/// Every `u_ii` is nonzero and not NaN, checked by `debug_assert!`; release behaviour as for [`solve_lower`].
#[inline]
pub fn solve_upper<S: Real, const N: usize>(u: &Matrix<S, N, N>, b: Vector<S, N>) -> Vector<S, N> {
    let mut x = [S::zero(); N];
    for i in (0..N).rev() {
        let uii = u.get(i, i);
        debug_assert!(uii.abs().value_f64() > 0.0, "solve_upper: zero diagonal");
        x[i] = (b.0[i] - sum((i + 1..N).map(|k| u.get(i, k) * x[k]))) / uii;
    }
    Vector(x)
}

/// Solves `l^T x = b` by back substitution, reading only the lower triangle of `l`, by column.
///
/// The operations of `solve_upper(&l.transpose(), b)` on the same operands, so the result is
/// bit-identical to it for every input, the sign and payload of a NaN from arithmetic aside
/// (`NUMERICS.md` §15.6, `PHASE2.md` §3): `x_i = (b_i - sum_{k>i} l_ki x_k) / l_ii`, the sum in
/// increasing `k`, `i` from `N - 1` down; `l_ki` for `k > i` is column `i` below its diagonal.
/// Domain and release behaviour as for [`solve_lower`].
#[inline]
pub(crate) fn solve_lower_transposed<S: Real, const N: usize>(
    l: &Matrix<S, N, N>,
    b: Vector<S, N>,
) -> Vector<S, N> {
    let mut x = [S::zero(); N];
    for i in (0..N).rev() {
        let lii = l.get(i, i);
        debug_assert!(
            lii.abs().value_f64() > 0.0,
            "solve_lower_transposed: zero diagonal"
        );
        x[i] = (b.0[i] - sum((i + 1..N).map(|k| l.get(k, i) * x[k]))) / lii;
    }
    Vector(x)
}

/// Solves `A x = b` from the Cholesky factor `l` of `A` (`L L^T x = b`): forward substitution
/// with [`solve_lower`], then back substitution against `l^T` read by column, with no
/// transposed copy.
///
/// Bit-identical to `solve_upper(&l.transpose(), solve_lower(&l, b))` for every input and every
/// `S`, the sign and payload of a NaN from arithmetic aside (`NUMERICS.md` §15.6, `PHASE2.md` §3);
/// that composition is the reference twin. Backward stable when the
/// factorization ran to completion: `(A + E) x = b`, `|E| <= gamma_{3N+1} |l| |l|^T` (Higham,
/// Thm 10.4).
///
/// # Domain
///
/// Every `l_ii` is nonzero and not NaN, checked by `debug_assert!` (as in [`solve_lower`], which
/// runs first). A release build never panics: a zero diagonal gives `inf` or NaN. The upper
/// triangle of `l` is not read. The factor of a `chol` whose mask was clear is a legal input, and
/// the result solves nothing.
///
/// # Example
///
/// ```
/// use helicoid_linalg::{chol, chol_solve, Mat2, Vector};
///
/// let a = Mat2::from_cols([Vector([4.0_f64, 2.0]), Vector([2.0, 5.0])]);
/// let (l, pd) = chol(&a);
/// // `a * (1, 2) = (8, 12)`.
/// let x = chol_solve(&l, Vector([8.0, 12.0]));
/// assert!(pd && x.0 == [1.0, 2.0]);
/// ```
#[inline]
pub fn chol_solve<S: Real, const N: usize>(l: &Matrix<S, N, N>, b: Vector<S, N>) -> Vector<S, N> {
    solve_lower_transposed(l, solve_lower(l, b))
}
