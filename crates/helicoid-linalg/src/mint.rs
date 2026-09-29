//! `From`/`Into` between the fixed-size types and `mint` (`docs/PHASE2.md` §7), for `f32` and
//! `f64` only.
//!
//! A conversion moves the scalars and never computes on them: bits, `-0` and NaN payloads
//! survive. Component `i` of a vector or point is mint's `x, y, z, w` in that order. A
//! [`Matrix`] is column-major, so mint's `ColumnMatrixN` field `x` is column 0, and entry
//! `(r, c)` is field `c`, component `r`. `mint` has no 4-D point (`Point2` and `Point3`
//! only), so `Point` converts at `N` = 2 and 3.

use crate::{Matrix, Point, Vector};

/// `From` both ways between `$ours<$s, $n>` and `$theirs<$s>` through the array `mint` already
/// converts to and from.
macro_rules! array_pair {
    ($ours:ident, $theirs:ident, $s:ty, $n:literal) => {
        impl From<::mint::$theirs<$s>> for $ours<$s, $n> {
            #[inline]
            fn from(v: ::mint::$theirs<$s>) -> Self {
                Self(<[$s; $n]>::from(v))
            }
        }

        impl From<$ours<$s, $n>> for ::mint::$theirs<$s> {
            #[inline]
            fn from(v: $ours<$s, $n>) -> Self {
                Self::from(v.0)
            }
        }
    };
}

macro_rules! matrix_pair {
    ($theirs:ident, $s:ty, $n:literal) => {
        impl From<::mint::$theirs<$s>> for Matrix<$s, $n, $n> {
            #[inline]
            fn from(m: ::mint::$theirs<$s>) -> Self {
                Self::from_cols(<[[$s; $n]; $n]>::from(m).map(Vector))
            }
        }

        impl From<Matrix<$s, $n, $n>> for ::mint::$theirs<$s> {
            #[inline]
            fn from(m: Matrix<$s, $n, $n>) -> Self {
                Self::from(::core::array::from_fn(|c| m.col(c).0))
            }
        }
    };
}

macro_rules! scalar {
    ($s:ty) => {
        array_pair!(Vector, Vector2, $s, 2);
        array_pair!(Vector, Vector3, $s, 3);
        array_pair!(Vector, Vector4, $s, 4);
        array_pair!(Point, Point2, $s, 2);
        array_pair!(Point, Point3, $s, 3);
        matrix_pair!(ColumnMatrix2, $s, 2);
        matrix_pair!(ColumnMatrix3, $s, 3);
        matrix_pair!(ColumnMatrix4, $s, 4);
    };
}

scalar!(f32);
scalar!(f64);
