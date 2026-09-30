//! Tests of the `mint` conversions: the component and column order against `mint`'s named
//! fields, then bitwise round trips both ways, so nothing is computed on the way. Specials
//! (`-0`, infinities, subnormals, signalling and payload NaNs) are pinned in every slot, since
//! random bit patterns never reach them. Values are compared with `to_bits`, never with `==`.

use crate::{Mat2, Mat3, Matrix, Point, Vector};
use proptest::prelude::*;

/// Distinct, exactly representable in `f32`, so the same literals serve both scalars.
const V: [f64; 16] = [
    1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
];

macro_rules! scalar_tests {
    ($m:ident, $s:ty, $bits:ty, $snan:expr, $qnan:expr, $sign:expr) => {
        mod $m {
            use super::*;

            fn v<const K: usize>() -> [$s; K] {
                core::array::from_fn(|i| V[i] as $s)
            }

            #[test]
            fn vector_component_order() {
                let a: Vector<$s, 4> = mint::Vector4 {
                    x: 1.0,
                    y: 2.0,
                    z: 3.0,
                    w: 4.0,
                }
                .into();
                assert_eq!(a.0.map(<$s>::to_bits), v::<4>().map(<$s>::to_bits));
                let m: mint::Vector3<$s> = Vector(v::<3>()).into();
                assert_eq!(
                    [m.x, m.y, m.z].map(<$s>::to_bits),
                    v::<3>().map(<$s>::to_bits)
                );
                let m: mint::Vector2<$s> = Vector(v::<2>()).into();
                assert_eq!([m.x.to_bits(), m.y.to_bits()], v::<2>().map(<$s>::to_bits));
            }

            #[test]
            fn point_component_order() {
                let p: Point<$s, 3> = mint::Point3 {
                    x: 1.0,
                    y: 2.0,
                    z: 3.0,
                }
                .into();
                assert_eq!(p.0.map(<$s>::to_bits), v::<3>().map(<$s>::to_bits));
                let m: mint::Point2<$s> = Point(v::<2>()).into();
                assert_eq!([m.x.to_bits(), m.y.to_bits()], v::<2>().map(<$s>::to_bits));
            }

            /// Entry `(r, c)` of ours is component `r` of mint's column `c` (`x, y, z, w`), and
            /// the conversion is not a transpose: the matrix is not symmetric.
            #[test]
            fn matrix_columns() {
                let m3 = mint::ColumnMatrix3 {
                    x: mint::Vector3 {
                        x: 1.0,
                        y: 2.0,
                        z: 3.0,
                    },
                    y: mint::Vector3 {
                        x: 4.0,
                        y: 5.0,
                        z: 6.0,
                    },
                    z: mint::Vector3 {
                        x: 7.0,
                        y: 8.0,
                        z: 9.0,
                    },
                };
                let a: Mat3<$s> = m3.into();
                for c in 0..3 {
                    for r in 0..3 {
                        assert_eq!(a.get(r, c).to_bits(), (V[3 * c + r] as $s).to_bits());
                    }
                }
                let back: mint::ColumnMatrix3<$s> = a.into();
                assert_eq!(
                    [back.x.y, back.y.x, back.z.z].map(<$s>::to_bits),
                    [2.0, 4.0, 9.0].map(<$s>::to_bits)
                );

                let a: Mat2<$s> = mint::ColumnMatrix2 {
                    x: mint::Vector2 { x: 1.0, y: 2.0 },
                    y: mint::Vector2 { x: 3.0, y: 4.0 },
                }
                .into();
                assert_eq!(a.get(1, 0).to_bits(), (2.0 as $s).to_bits());
                assert_eq!(a.get(0, 1).to_bits(), (3.0 as $s).to_bits());

                let a: Matrix<$s, 4, 4> =
                    mint::ColumnMatrix4::from(core::array::from_fn::<[$s; 4], 4, _>(|c| {
                        core::array::from_fn(|r| V[4 * c + r] as $s)
                    }))
                    .into();
                assert_eq!(a.get(3, 0).to_bits(), (4.0 as $s).to_bits());
                assert_eq!(a.get(0, 3).to_bits(), (13.0 as $s).to_bits());
            }

            /// Every scalar of `s` goes through every conversion; bits compared, never values.
            fn round_trip(s: [$s; 16]) -> Result<(), TestCaseError> {
                let b = |x: &[$s]| x.iter().map(|f| f.to_bits()).collect::<std::vec::Vec<_>>();

                let x = Vector([s[0], s[1]]);
                let y: Vector<$s, 2> = mint::Vector2::from(x).into();
                prop_assert_eq!(b(&x.0), b(&y.0));
                let x = Vector([s[0], s[1], s[2]]);
                let y: Vector<$s, 3> = mint::Vector3::from(x).into();
                prop_assert_eq!(b(&x.0), b(&y.0));
                let x = Vector([s[0], s[1], s[2], s[3]]);
                let y: Vector<$s, 4> = mint::Vector4::from(x).into();
                prop_assert_eq!(b(&x.0), b(&y.0));

                let x = Point([s[4], s[5]]);
                let y: Point<$s, 2> = mint::Point2::from(x).into();
                prop_assert_eq!(b(&x.0), b(&y.0));
                let x = Point([s[4], s[5], s[6]]);
                let y: Point<$s, 3> = mint::Point3::from(x).into();
                prop_assert_eq!(b(&x.0), b(&y.0));

                let cols = |n: usize| -> std::vec::Vec<$s> { s[..n * n].to_vec() };
                let m2 = Mat2::<$s>::from_cols(core::array::from_fn(|c| {
                    Vector(core::array::from_fn(|r| s[2 * c + r]))
                }));
                let r2: Mat2<$s> = mint::ColumnMatrix2::from(m2).into();
                prop_assert_eq!(b(&cols(2)), b(&flat(&r2)));
                let m3 = Mat3::<$s>::from_cols(core::array::from_fn(|c| {
                    Vector(core::array::from_fn(|r| s[3 * c + r]))
                }));
                let r3: Mat3<$s> = mint::ColumnMatrix3::from(m3).into();
                prop_assert_eq!(b(&cols(3)), b(&flat(&r3)));
                let m4 = Matrix::<$s, 4, 4>::from_cols(core::array::from_fn(|c| {
                    Vector(core::array::from_fn(|r| s[4 * c + r]))
                }));
                let r4: Matrix<$s, 4, 4> = mint::ColumnMatrix4::from(m4).into();
                prop_assert_eq!(b(&cols(4)), b(&flat(&r4)));
                Ok(())
            }

            /// Random bit patterns essentially never hit `-0`, an infinity or a specific NaN, so
            /// the specials are placed deterministically: rotating the pool puts each of them in
            /// each of the 16 slots.
            #[test]
            fn specials_in_every_slot() -> Result<(), TestCaseError> {
                let pool: [$s; 11] = [
                    0.0,
                    -0.0,
                    <$s>::INFINITY,
                    <$s>::NEG_INFINITY,
                    <$s>::MAX,
                    <$s>::MIN_POSITIVE,
                    <$s>::from_bits(1),
                    <$s>::from_bits(<$s>::MIN_POSITIVE.to_bits() - 1),
                    <$s>::from_bits($snan),
                    <$s>::from_bits($qnan),
                    <$s>::from_bits($qnan | $sign),
                ];
                for k in 0..pool.len() {
                    round_trip(core::array::from_fn(|i| pool[(i + k) % pool.len()]))?;
                }
                Ok(())
            }

            proptest! {
                #![proptest_config(ProptestConfig {
                    cases: 4096,
                    failure_persistence: None,
                    ..ProptestConfig::default()
                })]

                #[test]
                fn round_trips(bits in prop::collection::vec(any::<$bits>(), 16)) {
                    round_trip(core::array::from_fn(|i| <$s>::from_bits(bits[i])))?;
                }
            }
        }
    };
}

/// `m` in column-major order.
fn flat<S: crate::Real, const N: usize>(m: &Matrix<S, N, N>) -> std::vec::Vec<S> {
    (0..N * N).map(|i| m.get(i % N, i / N)).collect()
}

// A signalling NaN, a quiet NaN with a payload, and the sign bit.
scalar_tests!(
    f64_,
    f64,
    u64,
    0x7ff0_0000_0000_0001,
    0x7ff8_0000_0000_1234,
    1 << 63
);
scalar_tests!(f32_, f32, u32, 0x7f80_0001, 0x7fc0_1234, 1 << 31);
