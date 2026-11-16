//! Triangular 3x3 values (`docs/PHASE2.md` §4, decision `0117`): the factors `rq3` and `lq3`
//! return, holding the six entries their structure leaves.

use crate::matrix::Mat3;
use crate::real::{Blend, Real};
use crate::vector::{Vec3, Vector};

/// An upper triangular 3x3 matrix,
/// `[[u11, u12, u13], [0, u22, u23], [0, 0, u33]]`: `rq3`'s triangle.
///
/// The six entries are read and written in a named order (`docs/API.md` R3); the storage is not a
/// contract (D2).
///
/// ```
/// use helicoid_linalg::UpperTri3;
///
/// let u = UpperTri3::from_row_major([1.0_f64, 2.0, 3.0, 4.0, 5.0, 6.0]);
/// assert_eq!(u.to_matrix().get(1, 2), 5.0);
/// assert_eq!(u.to_matrix().get(2, 0), 0.0);
/// assert_eq!(u.transpose().to_row_major(), [1.0, 2.0, 4.0, 3.0, 5.0, 6.0]);
/// ```
#[derive(Clone, Copy, Debug)]
pub struct UpperTri3<S> {
    a: [S; 6],
}

/// A lower triangular 3x3 matrix,
/// `[[l11, 0, 0], [l21, l22, 0], [l31, l32, l33]]`: `lq3`'s triangle.
///
/// The six entries are read and written in a named order (`docs/API.md` R3); the storage is not a
/// contract (D2).
#[derive(Clone, Copy, Debug)]
pub struct LowerTri3<S> {
    a: [S; 6],
}

/// The transpose's permutation of the row-major entries, the same both ways:
/// `(u11, u12, u13, u22, u23, u33)` and `(l11, l21, l22, l31, l32, l33)`.
#[inline]
fn transposed<S: Copy>([a, b, c, d, e, f]: [S; 6]) -> [S; 6] {
    [a, b, d, c, e, f]
}

impl<S: Real> UpperTri3<S> {
    /// The identity.
    #[inline]
    pub fn identity() -> Self {
        let (z, o) = (S::zero(), S::one());
        Self {
            a: [o, z, z, o, z, o],
        }
    }

    /// From the entries row by row, `[u11, u12, u13, u22, u23, u33]`.
    #[inline]
    pub fn from_row_major(a: [S; 6]) -> Self {
        Self { a }
    }

    /// The entries row by row, `[u11, u12, u13, u22, u23, u33]`.
    #[inline]
    pub fn to_row_major(&self) -> [S; 6] {
        self.a
    }

    /// `(u11, u22, u33)`.
    #[inline]
    pub fn diagonal(&self) -> Vec3<S> {
        Vector([self.a[0], self.a[3], self.a[5]])
    }

    /// The matrix, with `+0` below the diagonal.
    #[inline]
    pub fn to_matrix(&self) -> Mat3<S> {
        let z = S::zero();
        let [a, b, c, d, e, f] = self.a;
        Mat3::from_rows([Vector([a, b, c]), Vector([z, d, e]), Vector([z, z, f])])
    }

    /// The transpose, a lower triangle; every entry keeps its bits.
    #[inline]
    pub fn transpose(&self) -> LowerTri3<S> {
        LowerTri3 {
            a: transposed(self.a),
        }
    }
}

impl<S: Real> LowerTri3<S> {
    /// The identity.
    #[inline]
    pub fn identity() -> Self {
        UpperTri3::identity().transpose()
    }

    /// From the entries row by row, `[l11, l21, l22, l31, l32, l33]`.
    #[inline]
    pub fn from_row_major(a: [S; 6]) -> Self {
        Self { a }
    }

    /// The entries row by row, `[l11, l21, l22, l31, l32, l33]`.
    #[inline]
    pub fn to_row_major(&self) -> [S; 6] {
        self.a
    }

    /// `(l11, l22, l33)`.
    #[inline]
    pub fn diagonal(&self) -> Vec3<S> {
        Vector([self.a[0], self.a[2], self.a[5]])
    }

    /// The matrix, with `+0` above the diagonal.
    #[inline]
    pub fn to_matrix(&self) -> Mat3<S> {
        let z = S::zero();
        let [a, b, c, d, e, f] = self.a;
        Mat3::from_rows([Vector([a, z, z]), Vector([b, c, z]), Vector([d, e, f])])
    }

    /// The transpose, an upper triangle; every entry keeps its bits.
    #[inline]
    pub fn transpose(&self) -> UpperTri3<S> {
        UpperTri3 {
            a: transposed(self.a),
        }
    }
}

impl<S: Real> Blend<S> for UpperTri3<S> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self {
            a: <[S; 6]>::blend(m, t.a, f.a),
        }
    }
}

impl<S: Real> Blend<S> for LowerTri3<S> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self {
            a: <[S; 6]>::blend(m, t.a, f.a),
        }
    }
}
