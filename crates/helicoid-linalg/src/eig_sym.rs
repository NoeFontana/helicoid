//! Eigendecomposition of a symmetric `N x N` matrix, `4 <= N <= 16` (`docs/NUMERICS.md` §20,
//! decision `0111`).
//!
//! `eig3`'s cyclic Jacobi at general `N`: `eig3`'s rotation, scaling window, deterministic exit,
//! polish, sort and Newton normalisation, with a relative rotation threshold. Every decision is a
//! mask, every square root and division has a safe argument, and the only operations are
//! `+ - * /`, `sqrt` and `abs`.

use crate::eig3::{angle, smallest_normal, unit, window, window_scale};
use crate::matrix::Matrix;
use crate::real::{is_finite, pow2, Blend, Mask, Precision, Real};
use crate::vector::Vector;

/// The most sweeps before the polish, `CAP[N - 4]` at `N = 4..=16`: one more than the most any
/// matrix of `0111`'s convergence study needed to leave (the sweep that rotates nothing counted),
/// at `f64`, which needs as many as `f32` or more at every `N`: 7 at `N = 4`, 11 at 9, 14 at 16, of
/// `2e4` matrices in each of 20 families and a hill climb on the entries.
pub(crate) const CAP: [usize; 13] = [8, 10, 11, 11, 12, 12, 12, 13, 13, 14, 14, 14, 15];

/// The symmetric working matrix by columns, `a[c][r]`, both triangles kept.
type Sym<S, const N: usize> = [[S; N]; N];

/// The state of the iteration: the matrix, the accumulated rotations by columns, and `r_i =
/// sqrt(|a_ii|)`, the scale of the relative threshold.
struct State<S: Real, const N: usize> {
    a: Sym<S, N>,
    v: Sym<S, N>,
    r: [S; N],
}

/// `angle`'s `(ta, c, s)` for `d = a_qq - a_pp` and `a_pq`, at every scale the relative threshold
/// reaches (`NUMERICS.md` §20). Where `|2 a_pq| >= sqrt(tiny)`, `w^2` and `d^2 + w^2` are normal
/// and this is `angle` itself; below, [`small`].
#[inline(always)]
fn rotation<S: Real>(d: S, apq: S) -> (S, S, S) {
    let w = S::lit(2.0) * apq;
    let lim = S::lit(match S::PRECISION {
        Precision::F64 => const { pow2(-511) },
        Precision::F32 => const { pow2(-63) },
    });
    let below = w.abs().lt(lim);
    // A scalar takes one arm, so the common one pays no latency for the test; a lane outside the
    // rare one passes it a zero, which scales without overflow.
    S::branch(
        below,
        move || small(d, S::select(below, apq, S::zero())),
        move || {
            let (ta, c, s, _, _) = angle(d, apq, w * w);
            (ta, c, s)
        },
    )
}

/// [`rotation`] where `|2 a_pq| < sqrt(tiny)`, out of line: rarely taken. `angle` reads `ta` from
/// `w^2`, which underflows here, and where `|d| < sqrt(tiny)` too it takes its safe arm and leaves
/// the pair unrotated: either loses a small eigenvalue the relative threshold keeps. So where `|d|
/// < sqrt(tiny)` the pair is scaled by `2^512` (`f64`), `2^64` (`f32`), exact, under which the
/// angle is invariant: then `|d|, |w| < 2`, and a rotating pair has `|a_pq| > tiny`, so `w^2 >
/// 2^-1018`, `2^-122`, and `angle` takes its exact arm. Where `w^2` is still below `tiny`, `ta = t
/// a_pq`, which squares nothing.
#[cold]
#[inline(never)]
fn small<S: Real>(d: S, apq: S) -> (S, S, S) {
    let (lim, up, down) = match S::PRECISION {
        Precision::F64 => (
            const { pow2(-511) },
            const { pow2(512) },
            const { pow2(-512) },
        ),
        Precision::F32 => (const { pow2(-63) }, const { pow2(64) }, const { pow2(-64) }),
    };
    let one = S::one();
    let near = d.abs().lt(S::lit(lim));
    let (up, down) = (
        S::select(near, S::lit(up), one),
        S::select(near, S::lit(down), one),
    );
    let (d, apq) = (d * up, apq * up);
    let w = S::lit(2.0) * apq;
    let w2 = w * w;
    let (ta, c, s, _, t) = angle(d, apq, w2);
    let ta = S::select(w2.lt(smallest_normal::<S>()), t * apq, ta);
    (ta * down, c, s)
}

impl<S: Real, const N: usize> State<S, N> {
    /// The rotation of `(p, q)` (`NUMERICS.md` §20) on columns and rows `p`, `q` of `a` and columns
    /// `p`, `q` of `v`. Where `keep` is set every entry is selected as it was, bit for bit: a
    /// rotation by `(1, 0)` would turn a `-0` beside a negative entry into `+0`. A scalar does not
    /// run a kept one at all. With `SCALE` the threshold's `r_p`, `r_q` follow; the polish reads
    /// neither.
    #[inline(always)]
    fn rotate<const SCALE: bool>(&mut self, p: usize, q: usize, keep: S::Mask) {
        // A lane mask runs every rotation with its lanes held; a scalar skips a kept one.
        if keep.decide(|| true, || false, |_, _, _| false) {
            return;
        }
        let (app, aqq, apq) = (self.a[p][p], self.a[q][q], self.a[p][q]);
        let (ta, c, s) = rotation(aqq - app, apq);
        let turn = |x: S, y: S| {
            (
                S::select(keep, x, c * x - s * y),
                S::select(keep, y, s * x + c * y),
            )
        };
        for k in 0..N {
            (self.a[p][k], self.a[q][k]) = turn(self.a[p][k], self.a[q][k]);
        }
        // Mirrored in a second pass: one pass that also stores the rows, strided, keeps LLVM from
        // vectorising the columns and measured 5-7% slower at `N = 9, 16`.
        for k in 0..N {
            self.a[k][p] = self.a[p][k];
            self.a[k][q] = self.a[q][k];
        }
        let left = S::select(keep, apq, S::zero());
        self.a[p][p] = S::select(keep, app, app - ta);
        self.a[q][q] = S::select(keep, aqq, aqq + ta);
        self.a[q][p] = left;
        self.a[p][q] = left;
        for k in 0..N {
            (self.v[p][k], self.v[q][k]) = turn(self.v[p][k], self.v[q][k]);
        }
        if SCALE {
            self.r[p] = self.a[p][p].abs().sqrt();
            self.r[q] = self.a[q][q].abs().sqrt();
        }
    }
}

/// The rotations of [`eig_sym`]: the diagonal left, the accumulated columns, the power of two that
/// undoes the scaling (`NaN` where an entry is not finite), and the sweeps run before the polish.
/// `REL` is the relative threshold; without it, the absolute `|a_pq| <= u |A|`, which only the
/// tests use, to show what the relative one buys.
#[inline(always)]
pub(crate) fn jacobi<S: Real, const N: usize, const REL: bool>(
    a: &Matrix<S, N, N>,
) -> ([S; N], Sym<S, N>, S, usize) {
    let zero = S::zero();
    // Only the lower triangle is read; a non-finite entry iterates on the zero matrix (`eig3`).
    let mut finite = zero.le(zero);
    let mut m = zero;
    for c in 0..N {
        for r in c..N {
            let x = a.get(r, c);
            finite = finite.and(is_finite(x));
            let x = x.abs();
            m = S::select(m.lt(x), x, m);
        }
    }
    let m = S::select(finite, m, zero);
    let one = S::one();
    // Scalars only cross the branch; a product with `1` is exact.
    let ((f, g), back) = window_scale(m, finite, window::<S>(), (one, one), |f, g| (f, g));
    let mut st = State {
        a: [[zero; N]; N],
        v: core::array::from_fn(|c| core::array::from_fn(|r| S::lit(f64::from(u8::from(r == c))))),
        r: [zero; N],
    };
    for c in 0..N {
        for r in c..N {
            let x = S::select(finite, a.get(r, c), zero) * f * g;
            st.a[c][r] = x;
            st.a[r][c] = x;
        }
        st.r[c] = st.a[c][c].abs().sqrt();
    }

    let u = S::lit(match S::PRECISION {
        Precision::F64 => const { pow2(-53) },
        Precision::F32 => const { pow2(-24) },
    });
    let tiny = smallest_normal::<S>();
    let absolute = u * (m * f * g);
    let cap = const { CAP[N - 4] };
    let (mut k, mut sweeps) = (0, 0);
    while k < cap {
        let mut calm = zero.le(zero);
        for p in 0..N {
            for q in p + 1..N {
                // Rotate where `|a_pq| > max(u r_p r_q, tiny)`: the floor keeps a pair of exactly
                // zero diagonals from rotating forever on an entry far below `u |A|`, and is the
                // bound below a rotating `a_pq` that `small` relies on.
                let t = if REL { u * st.r[p] * st.r[q] } else { absolute };
                let t = S::select(t.lt(tiny), tiny, t);
                let keep = st.a[p][q].abs().le(t);
                calm = calm.and(keep);
                st.rotate::<true>(p, q, keep);
            }
        }
        sweeps += 1;
        // A scalar leaves after a sweep that rotated nothing; lanes run on, each held above.
        k = calm.decide(|| cap, || k + 1, |_, _, next| next);
    }
    // The polish (`0084` decision 4 at `N`): every pair once more, rotated whatever its size.
    let never = zero.lt(zero);
    for p in 0..N {
        for q in p + 1..N {
            st.rotate::<false>(p, q, never);
        }
    }
    (core::array::from_fn(|k| st.a[k][k]), st.v, back, sweeps)
}

/// The eigenvalues, ascending, and an orthonormal basis of eigenvectors of the symmetric `N x N`
/// matrix `a`, `4 <= N <= 16`.
///
/// Returns `(values, vectors)`: `values.0[0] <= ... <= values.0[N - 1]`, column `k` of `vectors`
/// an eigenvector of `values.0[k]`, `vectors` orthonormal with `det = +1`. The sign of a column is
/// what the rotations produce; only the parity is fixed.
///
/// Cyclic Jacobi with a relative threshold and a deterministic exit (`docs/NUMERICS.md` §20,
/// `0111`): sweeps over the pairs `(p, q)`, `p < q` in lexicographic order, each rotation `eig3`'s
/// (`NUMERICS.md` §17), applied where `|a_pq| > u sqrt(|a_pp| |a_qq|)` (at least the smallest
/// normal number), a pair below the square root of the smallest normal scaled by an exact power of
/// two first, until a sweep applies none; then one sweep that rotates every pair, which moves
/// the values by at most the entries left and gives a `Dual` its first-order lanes. The values are
/// the diagonal left and the vectors the product of the rotations, sorted ascending by a stable
/// odd-even transposition network whose exchanges are counted into the first column's sign to keep
/// `det = +1`, each column normalised once by a Newton step. The matrix is first scaled by an exact
/// power of two into `eig3`'s window, and the values scaled back.
///
/// # Accuracy
///
/// Backward stable: the eigenvalues are within a few `u |A|` of the exact ones and column `k`
/// within a few `u |A| / gap_k` radians, `gap_k` its distance to the nearest other eigenvalue.
/// For a positive definite `A` it is also accurate relative to each eigenvalue (Demmel and
/// Veselić): with `D = diag(sqrt(a_ii))` and `kappa` the condition number of `D^-1 A D^-1`, the
/// eigenvalue `l_k` is within a few `u kappa l_k`, and its vector within a few `u kappa /
/// relgap_k`, `relgap_k = min_j |l_k - l_j| / sqrt(l_k l_j)`. A graded matrix, whose small
/// eigenvalues a tridiagonal QR solver returns only to `u |A|`, keeps its digits, down to a
/// diagonal spanning `1e-200` (`f64`) and `1e-30` (`f32`) in the corpus. The corpus ids
/// `eig_sym_n{4,9,16}` read both, per stratum (`0111`).
///
/// # Domain
///
/// Every input is legal and nothing is asserted (`docs/API.md` R6): only the lower triangle is
/// read, and the caller owns symmetry. A non-finite entry gives `NaN` eigenvalues and the identity
/// frame. `N` outside `4..=16` is a compile error; `N = 3` is [`eig3`](crate::eig3). Stable rustdoc
/// does not check the error code, so the example below, the same call at `N = 4`, is the positive
/// control:
///
/// ```compile_fail,E0080
/// let a = helicoid_linalg::Matrix::<f64, 3, 3>::identity();
/// let _ = helicoid_linalg::eig_sym(&a);
/// ```
///
/// ```compile_fail,E0080
/// let a = helicoid_linalg::Matrix::<f64, 17, 17>::identity();
/// let _ = helicoid_linalg::eig_sym(&a);
/// ```
///
/// A `Dual` result differentiates the rotations taken, the polish included: where the eigenvalues
/// are simple that is first-order perturbation theory, to the convergence of the exit. At a
/// repeated eigenvalue the lanes are finite but in general no directional derivative, as `eig3`'s.
///
/// # Example
///
/// ```
/// use helicoid_linalg::{eig_sym, Matrix, Vector};
///
/// // A tridiagonal 4x4 with 2 on the diagonal and -1 beside it: 2 - 2 cos(k pi / 5).
/// let a = Matrix::from_rows([
///     Vector([2.0_f64, -1.0, 0.0, 0.0]),
///     Vector([-1.0, 2.0, -1.0, 0.0]),
///     Vector([0.0, -1.0, 2.0, -1.0]),
///     Vector([0.0, 0.0, -1.0, 2.0]),
/// ]);
/// let (values, vectors) = eig_sym(&a);
/// for k in 0..4 {
///     let exact = 2.0 - 2.0 * (core::f64::consts::PI * (k + 1) as f64 / 5.0).cos();
///     assert!((values.0[k] - exact).abs() < 1e-14);
///     assert!((vectors.col(k).norm_sq() - 1.0).abs() < 1e-15);
/// }
/// ```
#[inline]
pub fn eig_sym<S: Real, const N: usize>(a: &Matrix<S, N, N>) -> (Vector<S, N>, Matrix<S, N, N>) {
    const { assert!(4 <= N && N <= 16, "eig_sym: N is 4..=16; N = 3 is eig3") };
    let (mut l, mut v, back, _) = jacobi::<S, N, true>(a);
    let mut sign = S::one();
    // Odd-even transposition: `N` passes of adjacent exchanges where the second is strictly less,
    // stable, so equal eigenvalues keep the order the rotations left.
    for pass in 0..N {
        let mut i = pass % 2;
        while i + 1 < N {
            let swap = l[i + 1].lt(l[i]);
            let (x, y) = (l[i], l[i + 1]);
            l[i] = S::select(swap, y, x);
            l[i + 1] = S::select(swap, x, y);
            let (vx, vy) = (v[i], v[i + 1]);
            v[i] = <[S; N]>::blend(swap, vy, vx);
            v[i + 1] = <[S; N]>::blend(swap, vx, vy);
            sign = S::select(swap, -sign, sign);
            i += 2;
        }
    }
    (
        Vector(l.map(|l| l * back)),
        Matrix::from_cols(core::array::from_fn(|k| {
            let col = Vector(v[k]);
            unit(if k == 0 { col.scale(sign) } else { col })
        })),
    )
}
