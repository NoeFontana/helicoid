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
    pub diag: Mat3<S>,
    /// `B₁, …, B_N`: `col[i]` is the dense block at block row `i + 1`, block column `0`.
    pub col: [Mat3<S>; N],
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

fn zero3<S: Real>() -> Mat3<S> {
    Matrix::from_cols([Vector([S::zero(); 3]); 3])
}

/// The `D x D` zero matrix, the scratch `sandwich` overwrites block by block.
fn zeros<S: Real, const D: usize>() -> Matrix<S, D, D> {
    Matrix::from_cols([Vector([S::zero(); D]); D])
}

/// The `3 x 3` block `(bi, bj)` of a dense matrix.
fn block<S: Real, const D: usize>(m: &Matrix<S, D, D>, bi: usize, bj: usize) -> Mat3<S> {
    Matrix::from_cols(array::from_fn(|c| {
        Vector(array::from_fn(|r| m.get(3 * bi + r, 3 * bj + c)))
    }))
}

fn put<S: Real, const D: usize>(m: &mut Matrix<S, D, D>, bi: usize, bj: usize, b: &Mat3<S>) {
    for c in 0..3 {
        for r in 0..3 {
            m.set(3 * bi + r, 3 * bj + c, b.get(r, c));
        }
    }
}

/// The same, into a strided view.
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
    /// The computed `det A` is nonzero and finite, checked by `debug_assert!` (a NaN fails it):
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
        // `0 < |det| < ∞`: false for `±0`, NaN and `±∞`. `|det|` is bound, not taken twice: for a
        // `Dual` scalar `abs` is a select and a negation per derivative lane.
        let m = det.abs();
        debug_assert!(
            S::zero().lt(m).and(m.lt(S::lit(f64::INFINITY))).all(),
            "SEn3Jac::inverse: det A is zero, NaN or infinite"
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
        let a = self.diag;
        let at = a.transpose();
        // The `N` transposes of `col`, taken once: the inner loop below would retake each of them
        // on every block row.
        let bt: [Mat3<S>; N] = self.col.map(|b| b.transpose());
        // `M = J Σ`: block row 0 is `A Σ_0l`, block row `i` is `B_i Σ_0l + A Σ_il`. The scratch
        // starts at zero, not at `Σ`: both loops write every one of the `(N + 1)²` blocks, so
        // seeding from a copy of their input would copy `D²` scalars that nothing reads.
        let mut m = zeros::<S, D>();
        for l in 0..=N {
            let s0 = block(cov, 0, l);
            put(&mut m, 0, l, &(a * s0));
            for (i, b) in self.col.iter().enumerate() {
                put(&mut m, i + 1, l, &(*b * s0 + a * block(cov, i + 1, l)));
            }
        }
        // `M Jᵀ`: block column 0 is `M_i0 Aᵀ`, block column `j` is `M_i0 B_jᵀ + M_ij Aᵀ`.
        let mut out = zeros();
        for i in 0..=N {
            let m0 = block(&m, i, 0);
            put(&mut out, i, 0, &(m0 * at));
            for (j, b) in bt.iter().enumerate() {
                put(&mut out, i, j + 1, &(m0 * *b + block(&m, i, j + 1) * at));
            }
        }
        out
    }
}
