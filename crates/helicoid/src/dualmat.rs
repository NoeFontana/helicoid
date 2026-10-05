//! The dual-matrix algebra of the SE_N(3) Jacobians (`docs/PHASE3.md` §5, `0005`,
//! `NUMERICS.md` §2.2).

use crate::sen3::SEn3Tangent;
use crate::traits::{Jac, Tangent};
use core::array;
use helicoid_linalg::{Blend, Mask, Mat3, Matrix, Real, StridedMut, Vector};

/// A Jacobian of SE_N(3), `A + εB`: the rotation-first dense matrix
///
/// ```text
/// [ A    0   …   0 ]
/// [ B₁   A         ]
/// [ ⋮        ⋱     ]
/// [ B_N          A ]
/// ```
///
/// stored as `diag = A` and `col = [B₁, …, B_N]` (`NUMERICS.md` §2.2). `Ad`, `ad`, `J_r`, `J_l`
/// and their inverses all have this shape, with `A` their SO(3) counterpart
/// (`docs/maths/se3.md` SE.11(e)). The shape is closed under [`mul`](Jac::mul) and
/// [`inverse`](Jac::inverse):
///
/// ```text
/// (A + εB)(C + εD) = AC + ε (B_i C + A D_i),      (A + εB)⁻¹ = A⁻¹ − ε A⁻¹ B_i A⁻¹
/// ```
///
/// A product costs `27 + 54N` multiplications, against `27 (N + 1)³` for the dense one (81 and
/// 216 for `N = 1`; 135 and 729 for `N = 2`), and the dense matrix has determinant `det(A)^(N+1)`
/// (SE.11). The transpose is block upper triangular and leaves the shape, so there is no
/// `transpose`: [`apply_transpose`](Jac::apply_transpose) and [`sandwich`](Jac::sandwich) are the
/// two operations that need it. A dense matrix is produced only by
/// [`write_dense`](Jac::write_dense) and `sandwich` (`docs/API.md` R5). There is no `PartialEq`,
/// `Add` or `Sub`.
#[derive(Clone, Copy, Debug)]
pub struct SEn3Jac<S, const N: usize> {
    /// `A`, the diagonal block, repeated `N + 1` times on the dense diagonal.
    ///
    /// Crate-private: a block is a representation of the dual matrix, not a value of it, so it is
    /// reached from outside through `write_dense` and `apply` (`0025` decision 5, `0028` option A).
    pub(crate) diag: Mat3<S>,
    /// `B₁, …, B_N`: `col[i]` is the dense block at block row `i + 1`, block column `0`.
    ///
    /// Crate-private, as [`Self::diag`] is.
    pub(crate) col: [Mat3<S>; N],
}

impl<S: Real, const N: usize> Blend<S> for SEn3Jac<S, N> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self {
            diag: Mat3::blend(m, t.diag, f.diag),
            col: <[Mat3<S>; N]>::blend(m, t.col, f.col),
        }
    }
}

/// `DOF` of the tangent, the dimension that `write_dense`, `sandwich` and the twins tie to.
pub(crate) const fn dof<S: Real, const N: usize>() -> usize {
    <SEn3Tangent<S, N> as Tangent<S>>::DOF
}

#[inline]
fn zero3<S: Real>() -> Mat3<S> {
    Matrix::from_cols([Vector([S::zero(); 3]); 3])
}

/// `J w` for the dense `J` of `A + εB` acting on a `D`-vector: block `0` is `A w_0`, block `i + 1`
/// is `B_i w_0 + A w_{i+1}`.
///
/// The same expression and the same grouping as [`apply`](Jac::apply) on a tangent, so the two agree
/// bit for bit; this one takes the flat vector a dense column or row already is, with no tangent to
/// read and write back.
#[inline]
fn act<S: Real, const N: usize, const D: usize>(
    a: &Mat3<S>,
    b: &[Mat3<S>; N],
    w: &Vector<S, D>,
) -> Vector<S, D> {
    // Entry by entry, holding no block vector and gathering nothing: `w` is read in place and each
    // output entry is finished before the next starts. The three products are summed left to right,
    // the association `Matrix`'s `Mul` and `Mul<Vector>` both use, and `B_i w_0` is added to
    // `A w_{i+1}` in the order `apply` adds them -- so every entry is the same products in the same
    // order as the block form, bit for bit.
    //
    // Measured against two alternatives at `D = 6, 9, 12` over `f64`, `f32` and `Dual<f64, 3>`,
    // release, x86_64, minimum of 7 runs of 200 000 calls: the block form through a `D x D` scratch
    // costs 545 ns (`f64`, `D = 6`), 2406 ns (`f64`, `D = 12`) and 8078 ns (`Dual`, `D = 12`); this
    // one costs 48, 448 and 2328. A middle form that built each block's `Vector` and gathered it was
    // faster still on `f64`/`f32` (365 and 327 ns at `D = 12`) but 17% to 20% *slower* on `Dual`,
    // which `0006`'s per-precision no-regress bar refuses.
    let dot = |m: &Mat3<S>, r: usize, off: usize| {
        m.get(r, 0) * w.0[off] + m.get(r, 1) * w.0[off + 1] + m.get(r, 2) * w.0[off + 2]
    };
    Vector(array::from_fn(|i| {
        let (k, r) = (i / 3, i % 3);
        if k == 0 {
            dot(a, r, 0)
        } else {
            dot(&b[k - 1], r, 0) + dot(a, r, 3 * k)
        }
    }))
}

/// The same, into a strided view.
#[inline]
fn put_view<S: Real>(out: &mut StridedMut<'_, S>, bi: usize, bj: usize, b: &Mat3<S>) {
    for c in 0..3 {
        for r in 0..3 {
            out.set(3 * bi + r, 3 * bj + c, b.get(r, c));
        }
    }
}

impl<S: Real, const N: usize> Jac<S, SEn3Tangent<S, N>> for SEn3Jac<S, N> {
    #[inline]
    fn identity() -> Self {
        Self {
            diag: Mat3::identity(),
            col: [zero3(); N],
        }
    }

    /// `(A + εB)(C + εD) = AC + ε (B_i C + A D_i)`: `27 + 54N` multiplications, and every entry
    /// a left-to-right sum of three products (or of two such sums, for a block of `col`).
    #[inline]
    fn mul(&self, o: &Self) -> Self {
        Self {
            diag: self.diag * o.diag,
            col: array::from_fn(|i| self.col[i] * o.diag + self.diag * o.col[i]),
        }
    }

    /// `A⁻¹ − ε A⁻¹ B_i A⁻¹`, `A⁻¹` from [`Mat3::inverse_adj`]: the inverse of the dense matrix
    /// (`docs/maths/se3.md` SE.11(b)), `54N` multiplications beyond `A⁻¹`.
    ///
    /// # Domain
    ///
    /// The computed `det A` is neither zero nor infinite, checked by `debug_assert!`:
    /// nonsingular to working precision, which is less than exactly invertible. The dense matrix
    /// has determinant `det(A)^(N+1)`. `inverse_adj` leaves what its `det` means to the caller,
    /// and this method reads it only for that assertion: a caller who must decide what a small
    /// `det A` means, for `J_r` near `θ = 2π` say, reads `self.diag.inverse_adj().1` first. A
    /// release build divides regardless: non-finite entries at `det A = 0`, and a finite, wrong
    /// zero where `det A` overflows. `adj` and `det` are not scaled, so the entries of `A` must
    /// lie in the range `inverse_adj` states (about `1e-100` to `1e100` for `f64`, `1e-12` to
    /// `1e12` for `f32`) or the quotient is lost: the assertion catches an underflowed and an
    /// overflowed `det`, not the digits lost short of either.
    #[inline]
    fn inverse(&self) -> Self {
        let (inv, det) = self.diag.inverse_adj();
        // "No lane has `det = 0`, and none has `|det| = ∞`" -- not "every lane has
        // `0 < |det| < ∞`", which the two spellings differ from only on NaN. A NaN determinant is
        // a NaN matrix, whose inverse is NaN: the answer a value function owes, not a panic
        // (`Mat3`'s `Jac::inverse` says the same, and `coeffs`'s `nonnegative` assert reads the
        // same way). An infinite one still fires, because `adj/det` is then a finite, wrong zero
        // and that is a domain violation with a finite input to blame. `|det|` is bound, not taken
        // twice: for a `Dual` scalar `abs` is a select and a negation per derivative lane.
        let m = det.abs();
        debug_assert!(
            !m.le(S::zero()).or(S::lit(f64::INFINITY).le(m)).any(),
            "SEn3Jac::inverse: det A is zero or infinite"
        );
        Self {
            diag: inv,
            col: self.col.map(|b| -(inv * b * inv)),
        }
    }

    #[inline]
    fn neg(&self) -> Self {
        Self {
            diag: -self.diag,
            col: self.col.map(|b| -b),
        }
    }

    /// `[A φ; B_i φ + A ρ_i]`, `9 + 18N` multiplications.
    #[inline]
    fn apply(&self, t: &SEn3Tangent<S, N>) -> SEn3Tangent<S, N> {
        SEn3Tangent {
            phi: self.diag * t.phi,
            rho: array::from_fn(|i| self.col[i] * t.phi + self.diag * t.rho[i]),
        }
    }

    /// `[Aᵀ φ + Σ_i B_iᵀ ρ_i; Aᵀ ρ_i]`, `9 + 18N` multiplications, the sum over `i` left to right
    /// after `Aᵀ φ`.
    #[inline]
    fn apply_transpose(&self, t: &SEn3Tangent<S, N>) -> SEn3Tangent<S, N> {
        let at = self.diag.transpose();
        let phi = self
            .col
            .iter()
            .zip(&t.rho)
            .fold(at * t.phi, |acc, (b, r)| acc + b.transpose() * *r);
        SEn3Tangent {
            phi,
            rho: t.rho.map(|r| at * r),
        }
    }

    /// Writes the `(3 + 3N)`-square matrix of the type's documentation, structural zeros
    /// included.
    ///
    /// # Domain
    ///
    /// `out` is `DOF x DOF`, checked by `debug_assert!`. A smaller view panics in the strided
    /// access, the one documented panic class (D11).
    #[inline]
    fn write_dense(&self, out: &mut StridedMut<'_, S>) {
        debug_assert!(
            out.rows() == dof::<S, N>() && out.cols() == dof::<S, N>(),
            "Jac::write_dense: the view is not DOF x DOF"
        );
        // The structure is fixed, so each block is written once from a reference: no `Option` of a
        // `Mat3` to construct and copy per entry. `bi == bj` is tested first, so `bj == 0` is
        // reached only with `bi > 0`.
        let zero = zero3();
        for bi in 0..=N {
            for bj in 0..=N {
                let blk = if bi == bj {
                    &self.diag
                } else if bj == 0 {
                    &self.col[bi - 1]
                } else {
                    &zero
                };
                put_view(out, bi, bj, blk);
            }
        }
    }

    /// `J Σ Jᵀ`, formed block by block from the structure: first `J Σ`, then `(J Σ) Jᵀ`, each a
    /// left-to-right sum of block products. `Σ` need not be symmetric, and the result is not
    /// symmetrized.
    ///
    /// `D` must be `3 + 3N`, asserted at monomorphization (`cargo build` and `cargo test`, not
    /// `cargo check`):
    ///
    /// ```compile_fail,E0080
    /// use helicoid::{Jac, SEn3Jac};
    /// use helicoid_linalg::Matrix;
    /// let j: SEn3Jac<f64, 1> = Jac::identity();
    /// let _ = j.sandwich::<7>(&Matrix::<f64, 7, 7>::identity());
    /// ```
    ///
    /// Positive control, so the failure above comes from the assertion and not from an import:
    ///
    /// ```
    /// use helicoid::{Jac, SEn3Jac};
    /// use helicoid_linalg::Matrix;
    /// let j: SEn3Jac<f64, 1> = Jac::identity();
    /// let _ = j.sandwich::<6>(&Matrix::<f64, 6, 6>::identity());
    /// ```
    #[inline]
    fn sandwich<const D: usize>(&self, cov: &Matrix<S, D, D>) -> Matrix<S, D, D> {
        const { assert!(D == dof::<S, N>()) };
        // `M = J Σ` by columns, then `M Jᵀ` by rows: row `r` of `M Jᵀ` is `J` applied to row `r` of
        // `M`, since `(M Jᵀ)[r][k] = Σ_q J[k][q] M[r][q]`. So neither `Aᵀ` nor the `N` transposes of
        // `col` are formed, and neither half needs a `D x D` scratch -- `from_cols`/`from_rows`
        // build each result directly, where the block form zeroed two of them (two `memset`s of
        // `D²` scalars in the release asm, both then fully overwritten) and reached every entry
        // through `Matrix::get`/`set`.
        //
        // Every entry is the same products summed in the same order as the block form -- `act` says
        // why -- so the values are bit-identical, D16 holds and the recorded bounds are unmoved
        // rather than re-measured. Checked directly over 150 000 random `(J, Σ)` at `f64`, `f32` and
        // `Dual<f64, 3>` for `N = 1, 2, 3`: not one differing entry.
        let m = Matrix::from_cols(array::from_fn(|c| {
            act::<S, N, D>(&self.diag, &self.col, &cov.col(c))
        }));
        Matrix::from_rows(array::from_fn(|r| {
            act::<S, N, D>(&self.diag, &self.col, &m.row(r))
        }))
    }
}
