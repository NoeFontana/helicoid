//! Fixed-size column-major matrices (`docs/PHASE2.md` §4).
//!
//! Every product entry is `sum_k a_rk b_kc`, added in `k` order from the first term
//! ([`crate::vector`]'s reduction), for matrix-matrix and matrix-vector alike.

use crate::real::{Blend, Real};
use crate::vector::{sum, Vector};
use core::array;
use core::ops::{Add, Mul, Neg, Sub};

/// An `R x C` matrix, column-major: entry `(r, c)` is `cols[c][r]`. Layout is not a contract (D2).
///
/// There is no `PartialEq`, for the reason given on [`Vector`]. This does not compile:
///
/// ```compile_fail,E0369
/// use helicoid_linalg::Mat3;
/// fn same(a: Mat3<f64>, b: Mat3<f64>) -> bool { a == b }
/// ```
///
/// Positive control:
///
/// ```
/// use helicoid_linalg::Mat3;
/// fn prod(a: Mat3<f64>, b: Mat3<f64>) -> Mat3<f64> { a * b }
/// assert_eq!(prod(Mat3::identity(), Mat3::identity()).get(1, 1).to_bits(), 1.0_f64.to_bits());
/// ```
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct Matrix<S, const R: usize, const C: usize> {
    cols: [[S; R]; C],
}

/// A 2x2 matrix.
pub type Mat2<S> = Matrix<S, 2, 2>;
/// A 3x3 matrix.
pub type Mat3<S> = Matrix<S, 3, 3>;

impl<S: Real, const R: usize, const C: usize> Matrix<S, R, C> {
    /// The matrix whose `c`-th column is `cols[c]`.
    #[inline]
    pub fn from_cols(cols: [Vector<S, R>; C]) -> Self {
        Self {
            cols: cols.map(|v| v.0),
        }
    }

    /// The matrix whose `r`-th row is `rows[r]`.
    #[inline]
    pub fn from_rows(rows: [Vector<S, C>; R]) -> Self {
        Self {
            cols: array::from_fn(|c| array::from_fn(|r| rows[r].0[c])),
        }
    }

    /// Entry `(r, c)`.
    ///
    /// # Domain
    ///
    /// `r < R` and `c < C`, checked by `debug_assert!`. A release build never panics and returns an
    /// unspecified value (in practice `+0`) outside.
    #[inline]
    pub fn get(&self, r: usize, c: usize) -> S {
        debug_assert!(r < R && c < C, "Matrix::get: index out of range");
        self.cols
            .get(c)
            .and_then(|col| col.get(r))
            .copied()
            .unwrap_or_else(S::zero)
    }

    /// Set entry `(r, c)` to `v`.
    ///
    /// # Domain
    ///
    /// `r < R` and `c < C`, checked by `debug_assert!`. A release build never panics and leaves the matrix
    /// unspecified (in practice unchanged) when outside.
    #[inline]
    pub fn set(&mut self, r: usize, c: usize, v: S) {
        debug_assert!(r < R && c < C, "Matrix::set: index out of range");
        if let Some(e) = self.cols.get_mut(c).and_then(|col| col.get_mut(r)) {
            *e = v;
        }
    }

    /// Column `c`.
    ///
    /// # Domain
    ///
    /// `c < C`, checked by `debug_assert!`. A release build never panics and returns an unspecified
    /// vector (in practice zero) outside.
    #[inline]
    pub fn col(&self, c: usize) -> Vector<S, R> {
        debug_assert!(c < C, "Matrix::col: index out of range");
        Vector(array::from_fn(|r| self.get(r, c)))
    }

    /// Row `r`.
    ///
    /// # Domain
    ///
    /// `r < R`, checked by `debug_assert!`. A release build never panics and returns an unspecified
    /// vector (in practice zero) outside.
    #[inline]
    pub fn row(&self, r: usize) -> Vector<S, C> {
        debug_assert!(r < R, "Matrix::row: index out of range");
        Vector(array::from_fn(|c| self.get(r, c)))
    }

    /// The transpose.
    #[inline]
    pub fn transpose(&self) -> Matrix<S, C, R> {
        Matrix {
            cols: array::from_fn(|r| array::from_fn(|c| self.cols[c][r])),
        }
    }

    /// `s * self`, entrywise.
    #[inline]
    pub fn scale(self, s: S) -> Self {
        Self {
            cols: self.cols.map(|col| col.map(|x| x * s)),
        }
    }
}

impl<S: Real, const N: usize> Matrix<S, N, N> {
    /// The identity: `+1` on the diagonal, `+0` elsewhere.
    #[inline]
    pub fn identity() -> Self {
        Self {
            cols: array::from_fn(|c| array::from_fn(|r| if r == c { S::one() } else { S::zero() })),
        }
    }
}

impl<S: Real> Mat3<S> {
    /// `(adj(self) / det, det)`, with `det = c0 . (c1 x c2)` and `adj` the transpose of the
    /// cofactor matrix, whose rows are `c1 x c2`, `c2 x c0`, `c0 x c1` (`c_i` the columns), so
    /// `self * adj = det I`. Each entry is one division by `det`, no reciprocal.
    ///
    /// # Domain
    ///
    /// The caller decides what `det` means: for `det = 0` the entries are non-finite, and nothing
    /// is asserted, since a singular matrix is legal input.
    ///
    /// `adj` and `det` are quadratic and cubic in the entries and are not scaled, so they must be
    /// normal finite numbers. For a well-conditioned matrix of entries of magnitude `m` that is
    /// about `1e-100 < m < 1e100` (`f64`) and `1e-12 < m < 1e12` (`f32`); beyond it `det`
    /// overflows to `inf` and the quotient is a finite, wrong `0`, or `det` underflows to `0`
    /// and the quotient is `inf`, or digits are lost.
    #[inline]
    pub fn inverse_adj(&self) -> (Self, S) {
        let (c0, c1, c2) = (self.col(0), self.col(1), self.col(2));
        let (r0, r1, r2) = (c1.cross(c2), c2.cross(c0), c0.cross(c1));
        let det = c0.dot(r0);
        let row = |r: Vector<S, 3>| Vector(r.0.map(|x| x / det));
        (Self::from_rows([row(r0), row(r1), row(r2)]), det)
    }
}

impl<S: Real, const R: usize, const C: usize> Add for Matrix<S, R, C> {
    type Output = Self;
    #[inline]
    fn add(self, o: Self) -> Self {
        Self {
            cols: array::from_fn(|c| array::from_fn(|r| self.cols[c][r] + o.cols[c][r])),
        }
    }
}

impl<S: Real, const R: usize, const C: usize> Sub for Matrix<S, R, C> {
    type Output = Self;
    #[inline]
    fn sub(self, o: Self) -> Self {
        Self {
            cols: array::from_fn(|c| array::from_fn(|r| self.cols[c][r] - o.cols[c][r])),
        }
    }
}

impl<S: Real, const R: usize, const C: usize> Neg for Matrix<S, R, C> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self {
            cols: self.cols.map(|col| col.map(Neg::neg)),
        }
    }
}

impl<S: Real, const R: usize, const K: usize, const C: usize> Mul<Matrix<S, K, C>>
    for Matrix<S, R, K>
{
    type Output = Matrix<S, R, C>;
    #[inline]
    fn mul(self, o: Matrix<S, K, C>) -> Matrix<S, R, C> {
        Matrix {
            cols: array::from_fn(|c| {
                array::from_fn(|r| sum((0..K).map(|k| self.cols[k][r] * o.cols[c][k])))
            }),
        }
    }
}

impl<S: Real, const R: usize, const C: usize> Mul<Vector<S, C>> for Matrix<S, R, C> {
    type Output = Vector<S, R>;
    #[inline]
    fn mul(self, v: Vector<S, C>) -> Vector<S, R> {
        Vector(array::from_fn(|r| {
            sum((0..C).map(|k| self.cols[k][r] * v.0[k]))
        }))
    }
}

impl<S: Real, const R: usize, const C: usize> Blend<S> for Matrix<S, R, C> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self {
            cols: <[[S; R]; C]>::blend(m, t.cols, f.cols),
        }
    }
}
