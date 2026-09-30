//! Tests of `chol`, `solve_lower`, `solve_upper`, `chol_solve`. Small cases are exact (integers
//! whose factor is an integer); the algebra runs under proptest for `f64`, `f32` and
//! `Dual<f64, 2>` on the pool of [`crate::linalg_tests`], with its shadow method: an identity holds
//! entrywise, per lane, to `bound(lane) * scale`, `scale` being the same computation on the
//! absolute values.
//!
//! The rounding bounds are Higham's (*Accuracy and Stability of Numerical Algorithms*): Thm 10.3
//! (`L L^T = A + E`, `|E| <= gamma_{N+1} |L||L^T|`), Thm 8.5 (a triangular solve is `(T + E) x =
//! b`, `|E| <= gamma_N |T|`) and Thm 10.4 (`(A + E) x = b`, `|E| <= gamma_{3N+1} |L||L^T|`), each
//! plus the roundings of the residual's own product. A derivative lane adds `extra` roundings per
//! level, measured. A test's doc records the worst `|diff| / (bound * scale)` for `f64`, `f32` and
//! `D`, all lanes, over 10^6 cases from a fixed seed:
//!
//! ```text
//! PROPTEST_CASES=1000000 PROPTEST_RNG_SEED=1 cargo nextest run --release -p helicoid-linalg chol_tests
//! ```

use crate::chol::solve_lower_transposed;
use crate::linalg_tests::{cfg, extra, flat, gamma, pool, ratio, shadow_m, shadow_v, unit, within};
use crate::linalg_tests::{Lane, Pool, D};
use crate::{chol, chol_solve, solve_lower, solve_upper, Mask, Mat3, Matrix, Vector};
use core::array;
use proptest::prelude::*;

/// The matrix with rows `rows`, every entry a constant of `S`.
fn mat<S: Lane, const N: usize>(rows: [[f64; N]; N]) -> Matrix<S, N, N> {
    Matrix::from_rows(rows.map(|r| Vector(r.map(|v| S::make(v, [0.0; 2])))))
}

fn vecn<S: Lane, const N: usize>(v: [f64; N]) -> Vector<S, N> {
    Vector(v.map(|v| S::make(v, [0.0; 2])))
}

fn finite<S: Lane>(x: &[S]) -> bool {
    x.iter()
        .all(|s| (0..S::LANES).all(|k| s.lane(k).is_finite()))
}

/// `|L||L^T|` in `D` arithmetic: the shadow scale of `L L^T`.
fn llt_scale<S: Lane, const N: usize>(l: &Matrix<S, N, N>) -> Matrix<D, N, N> {
    let s = shadow_m(l);
    s * s.transpose()
}

/// `L L^T - A` against `gamma_{2N+1} |L||L^T|`: Thm 10.3's `gamma_{N+1}` for the factor and
/// `gamma_N` for the product that reads it back. Only meaningful when the mask is set.
fn reconstructs<S: Lane, const N: usize>(
    a: &Matrix<S, N, N>,
    l: &Matrix<S, N, N>,
) -> Result<(), TestCaseError> {
    let bound = |k| gamma::<S>(2 * N + 1 + extra(k, 4 * N));
    within(
        "L L^T",
        ratio(
            &flat(&(*l * l.transpose())),
            &flat(a),
            &flat(&llt_scale(l)),
            &bound,
        ),
    )
}

fn strictly_upper_is_plus_zero<S: Lane, const N: usize>(l: &Matrix<S, N, N>) -> bool {
    (0..N).all(|c| (0..c).all(|r| (0..S::LANES).all(|k| l.get(r, c).lane(k).to_bits() == 0)))
}

/// `B B^T + I` from the pool: positive definite with `lambda_min >= 1`, so the mask must be set.
/// Measured worst: 0.80 (f64, D) and 0.80 (f32), at `N = 2`; the ratio falls with `N`.
fn spd<S: Lane, const N: usize>(pl: &[[f64; 3]]) -> Result<(), TestCaseError> {
    let b = Pool(pl.iter()).mat::<S, N, N>();
    let a = b * b.transpose() + Matrix::identity();
    let (l, pd) = chol(&a);
    prop_assert!(pd.all(), "SPD input reported as not positive definite");
    prop_assert!(strictly_upper_is_plus_zero(&l));
    reconstructs(&a, &l)
}

fn spd_all<S: Lane>(pl: &[[f64; 3]]) -> Result<(), TestCaseError> {
    spd::<S, 1>(pl)?;
    spd::<S, 2>(pl)?;
    spd::<S, 3>(pl)?;
    spd::<S, 4>(pl)?;
    spd::<S, 5>(pl)?;
    spd::<S, 6>(pl)
}

/// `M + M^T` from the pool, mostly indefinite: the output is finite whatever the mask, `diag(L)`
/// is positive, and a set mask means the factor reads back within the bound. `D` is checked for its value lane only,
/// which is bitwise the plain evaluation of `f64` (the mask included).
fn symmetric<S: Lane, const N: usize>(pl: &[[f64; 3]]) -> Result<(), TestCaseError> {
    let m = Pool(pl.iter()).mat::<S, N, N>();
    let a = m + m.transpose();
    let (l, pd) = chol(&a);
    prop_assert!(finite(&flat(&l)), "non-finite factor, mask {}", pd.all());
    prop_assert!(strictly_upper_is_plus_zero(&l));
    prop_assert!(
        (0..N).all(|j| l.get(j, j).lane(0) > 0.0),
        "L_jj is a positive root or 1"
    );
    if pd.all() {
        reconstructs(&a, &l)?;
    }
    Ok(())
}

/// Symmetric matrices of mantissa `m` in `[-1, 1]` times `2^e`, `e` over the whole finite range
/// (`f32` and `D` use `e / scale` so the entries stay finite in the scalar): squares and products
/// overflow routinely, and the factor must still be finite with a positive diagonal.
fn wide<S: Lane, const N: usize>(v: &[(f64, i32)], scale: i32) -> Result<(), TestCaseError> {
    let entry = |r: usize, c: usize| {
        let (m, e) = v[r.max(c) * N + r.min(c)];
        S::make(libm::ldexp(m, e / scale), [0.0; 2])
    };
    let a =
        Matrix::<S, N, N>::from_cols(array::from_fn(|c| Vector(array::from_fn(|r| entry(r, c)))));
    let (l, pd) = chol(&a);
    prop_assert!(finite(&flat(&l)), "non-finite factor, mask {}", pd.all());
    prop_assert!(strictly_upper_is_plus_zero(&l));
    prop_assert!((0..N).all(|j| l.get(j, j).lane(0) > 0.0));
    Ok(())
}

fn symmetric_all<S: Lane>(pl: &[[f64; 3]]) -> Result<(), TestCaseError> {
    symmetric::<S, 1>(pl)?;
    symmetric::<S, 2>(pl)?;
    symmetric::<S, 3>(pl)?;
    symmetric::<S, 4>(pl)?;
    symmetric::<S, 5>(pl)?;
    symmetric::<S, 6>(pl)
}

/// A dual matrix's value lane factors exactly like the plain matrix, mask included.
fn dual_value_is_plain<const N: usize>(pl: &[[f64; 3]]) -> Result<(), TestCaseError> {
    let m = Pool(pl.iter()).mat::<D, N, N>();
    let a = m + m.transpose();
    let plain = Matrix::<f64, N, N>::from_cols(array::from_fn(|c| {
        Vector(array::from_fn(|r| a.get(r, c).v))
    }));
    let ((ld, pd_d), (lp, pd_p)) = (chol(&a), chol(&plain));
    prop_assert!(pd_d == pd_p);
    for (x, y) in flat(&ld).iter().zip(flat(&lp)) {
        prop_assert!(x.v.to_bits() == y.to_bits());
    }
    Ok(())
}

/// A negative diagonal entry is never positive definite: the mask is clear and the factor finite.
fn negative_diagonal<S: Lane, const N: usize>(pl: &[[f64; 3]]) -> Result<(), TestCaseError> {
    let b = Pool(pl.iter()).mat::<S, N, N>();
    let spd = b * b.transpose() + Matrix::identity();
    for i in 0..N {
        let mut a = spd;
        a.set(i, i, -a.get(i, i));
        let (l, pd) = chol(&a);
        prop_assert!(!pd.any(), "negative a_{i}{i} reported positive definite");
        prop_assert!(finite(&flat(&l)));
    }
    Ok(())
}

/// `v v^T + eps I`, from exactly singular (`eps = 0`) to well posed. The mask may go either way
/// near the boundary; what must hold at every `eps` is a finite factor, and a set mask implies the
/// backward-stability bound.
fn near_singular<S: Lane>(pl: &[[f64; 3]]) -> Result<(), TestCaseError> {
    let v = Pool(pl.iter()).vec::<S, 4>();
    let vvt = Matrix::from_cols([v]) * Matrix::from_cols([v]).transpose();
    for e in (0..=64).step_by(4).map(Some).chain([None]) {
        let eps = e.map_or(0.0, |e| libm::ldexp(1.0, -e));
        let a = vvt + Matrix::identity().scale(S::make(eps, [0.0; 2]));
        let (l, pd) = chol(&a);
        prop_assert!(finite(&flat(&l)), "non-finite factor at eps = {eps}");
        if pd.all() {
            reconstructs(&a, &l)?;
        }
    }
    Ok(())
}

/// `(lower, upper)` residuals of `solve_lower` and `solve_upper` against `gamma_{2N} |T||x|`, with a
/// unit-scale-or-larger diagonal (`[1, 2)`) and `b` a rounded `T x`: Thm 8.5's `gamma_N` plus
/// `gamma_N` for the product that forms the residual. The upper solve reads the transpose.
/// Measured worst: 0.91 (f64, D) and 0.90 (f32), at `N = 2`.
fn triangular<S: Lane, const N: usize>(pl: &[[f64; 3]]) -> Result<(), TestCaseError> {
    let mut p = Pool(pl.iter());
    let raw = p.mat::<S, N, N>();
    let l = Matrix::from_cols(array::from_fn(|c| {
        Vector(array::from_fn(|r| {
            let e = raw.get(r, c);
            match r.cmp(&c) {
                core::cmp::Ordering::Less => S::zero(),
                core::cmp::Ordering::Equal => S::one() + e.abs(),
                core::cmp::Ordering::Greater => e,
            }
        }))
    }));
    let u = l.transpose();
    let x = p.vec::<S, N>();
    let bound = |k| gamma::<S>(2 * N + extra(k, 2 * N));
    let zero = [S::zero(); N];

    let xl = solve_lower(&l, l * x);
    let scale = shadow_m(&l) * shadow_v(xl);
    within("lower", ratio(&(l * xl - l * x).0, &zero, &scale.0, &bound))?;

    let xu = solve_upper(&u, u * x);
    let scale = shadow_m(&u) * shadow_v(xu);
    within("upper", ratio(&(u * xu - u * x).0, &zero, &scale.0, &bound))
}

fn triangular_all<S: Lane>(pl: &[[f64; 3]]) -> Result<(), TestCaseError> {
    triangular::<S, 1>(pl)?;
    triangular::<S, 2>(pl)?;
    triangular::<S, 3>(pl)?;
    triangular::<S, 4>(pl)?;
    triangular::<S, 5>(pl)?;
    triangular::<S, 6>(pl)
}

/// `A x = b` through `chol`, `solve_lower`, `solve_upper` on `B B^T + I`: the residual is within
/// `gamma_{4N+2} |L||L^T||x|`, Thm 10.4's `gamma_{3N+1}` plus `gamma_N` for `A x` and the
/// `gamma_{N+1}` gap between `|A|` and `|L||L^T|`. Measured worst: 0.67 (f64, D) and 0.67 (f32),
/// at `N = 1`.
fn solve_spd<S: Lane, const N: usize>(pl: &[[f64; 3]]) -> Result<(), TestCaseError> {
    let mut p = Pool(pl.iter());
    let b = p.mat::<S, N, N>();
    let a = b * b.transpose() + Matrix::identity();
    let (l, pd) = chol(&a);
    prop_assert!(pd.all());
    let rhs = a * p.vec::<S, N>();
    let x = chol_solve(&l, rhs);
    let bound = |k| gamma::<S>(4 * N + 2 + extra(k, 6 * N));
    let scale = llt_scale(&l) * shadow_v(x);
    within(
        "A x - b",
        ratio(&(a * x - rhs).0, &[S::zero(); N], &scale.0, &bound),
    )
}

fn solve_spd_all<S: Lane>(pl: &[[f64; 3]]) -> Result<(), TestCaseError> {
    solve_spd::<S, 1>(pl)?;
    solve_spd::<S, 2>(pl)?;
    solve_spd::<S, 3>(pl)?;
    solve_spd::<S, 4>(pl)?;
    solve_spd::<S, 5>(pl)?;
    solve_spd::<S, 6>(pl)
}

fn bits<S: Lane>(x: &[S]) -> std::vec::Vec<u64> {
    x.iter().map(|s| s.lane(0).to_bits()).collect()
}

/// Every lane's bits, not the value lane's only, with every NaN as one value: the sign and payload
/// of a NaN from arithmetic are unspecified in Rust and a release build may commute operands
/// (`NaN * -NaN`), as in `dual_value_is_plain_value` (`docs/PHASE2.md` §3). `+-0` and `+-inf` stay
/// distinct.
fn all_bits<S: Lane>(x: &[S]) -> std::vec::Vec<u64> {
    x.iter()
        .flat_map(|s| {
            (0..S::LANES).map(move |k| {
                let v = s.lane(k);
                if v.is_nan() {
                    f64::NAN.to_bits()
                } else {
                    v.to_bits()
                }
            })
        })
        .collect()
}

/// A raw cell: a selector per lane and a plain value. Selectors below 9 pick a special (`+-0`,
/// `+-inf`, `+-NaN`, the extremes of the range, a subnormal); the rest keep the plain value, so
/// every rounding, overflow and NaN path is reachable.
type Cell = ([u8; 3], f64);

fn special(sel: u8, v: f64) -> f64 {
    match sel {
        0 => 0.0,
        1 => -0.0,
        2 => f64::INFINITY,
        3 => f64::NEG_INFINITY,
        4 => f64::NAN,
        5 => -f64::NAN,
        6 => f64::MAX,
        7 => f64::MIN_POSITIVE,
        8 => f64::from_bits(1),
        _ => v,
    }
}

fn cell<S: Lane>((sel, v): Cell) -> S {
    S::make(
        special(sel[0], v),
        [special(sel[1], 2.0 * v), special(sel[2], -3.0 * v)],
    )
}

/// `n` cells whose selectors are drawn from `sel`: `0..24` for a wild mix, `9..24` for plain
/// values only.
fn cells(n: usize, sel: core::ops::Range<u8>) -> impl Strategy<Value = std::vec::Vec<Cell>> {
    prop::collection::vec((prop::array::uniform3(sel), -1e3f64..1e3), n)
}

/// A nonzero, non-NaN diagonal entry: the domain of the solves. A cell outside it becomes `2`.
fn diagonal<S: Lane>(c: Cell) -> S {
    let x = cell::<S>(c);
    if x.abs().value_f64() > 0.0 {
        x
    } else {
        S::make(2.0, [c.1, 0.0])
    }
}

/// `chol_solve` and the transposed solve against `solve_upper(&l.transpose(), ..)`, all lanes,
/// bit for bit (a NaN is one value, see `all_bits`). The lower triangle and `b` take every
/// special; the diagonal is in the domain; the strictly upper triangle is junk (never read).
fn solve_bits<S: Lane, const N: usize>(cells: &[Cell], b: &[Cell]) -> Result<(), TestCaseError> {
    let l = Matrix::<S, N, N>::from_cols(array::from_fn(|c| {
        Vector(array::from_fn(|r| match r.cmp(&c) {
            core::cmp::Ordering::Equal => diagonal(cells[r * N + c]),
            _ => cell(cells[r * N + c]),
        }))
    }));
    let b = Vector::<S, N>(array::from_fn(|i| cell(b[i])));
    let u = l.transpose();
    let y = solve_lower(&l, b);
    prop_assert_eq!(
        all_bits(&solve_lower_transposed(&l, y).0),
        all_bits(&solve_upper(&u, y).0)
    );
    prop_assert_eq!(
        all_bits(&chol_solve(&l, b).0),
        all_bits(&solve_upper(&u, solve_lower(&l, b)).0)
    );
    Ok(())
}

fn solve_bits_all<S: Lane>(cells: &[Cell], b: &[Cell]) -> Result<(), TestCaseError> {
    solve_bits::<S, 1>(cells, b)?;
    solve_bits::<S, 2>(cells, b)?;
    solve_bits::<S, 3>(cells, b)?;
    solve_bits::<S, 4>(cells, b)?;
    solve_bits::<S, 5>(cells, b)?;
    solve_bits::<S, 6>(cells, b)
}

fn known_factor<S: Lane>() {
    // The classic integer example: `L = [[2, 0, 0], [6, 1, 0], [-8, 5, 3]]`, `x = (1, 2, 3)`.
    let a = mat::<S, 3>([
        [4.0, 12.0, -16.0],
        [12.0, 37.0, -43.0],
        [-16.0, -43.0, 98.0],
    ]);
    let (l, pd) = chol(&a);
    let want = mat::<S, 3>([[2.0, 0.0, 0.0], [6.0, 1.0, 0.0], [-8.0, 5.0, 3.0]]);
    assert!(pd.all());
    assert_eq!(bits(&flat(&l)), bits(&flat(&want)));
    let x = vecn::<S, 3>([1.0, 2.0, 3.0]);
    let rhs = a * x;
    let y = solve_lower(&l, rhs);
    assert_eq!(bits(&y.0), bits(&(l.transpose() * x).0));
    assert_eq!(bits(&solve_upper(&l.transpose(), y).0), bits(&x.0));
    assert_eq!(bits(&solve_lower(&l, l * x).0), bits(&x.0));
    assert_eq!(bits(&chol_solve(&l, rhs).0), bits(&x.0));
}

#[test]
fn factors_a_known_integer_matrix_exactly() {
    known_factor::<f64>();
    known_factor::<f32>();
    known_factor::<D>();
}

#[test]
fn reads_only_the_lower_triangle() {
    let sym = mat::<f64, 3>([
        [4.0, 12.0, -16.0],
        [12.0, 37.0, -43.0],
        [-16.0, -43.0, 98.0],
    ]);
    let junk = mat::<f64, 3>([
        [4.0, f64::NAN, 1e300],
        [12.0, 37.0, f64::INFINITY],
        [-16.0, -43.0, 98.0],
    ]);
    let (a, b) = (chol(&sym), chol(&junk));
    assert!(a.1 && b.1);
    assert_eq!(bits(&flat(&a.0)), bits(&flat(&b.0)));
    // The solves read one triangle each.
    let t = mat::<f64, 2>([[2.0, f64::NAN], [3.0, 4.0]]);
    assert_eq!(
        bits(&solve_lower(&t, vecn([2.0, 7.0])).0),
        bits(&[1.0, 1.0])
    );
    let t = mat::<f64, 2>([[2.0, 3.0], [f64::NAN, 4.0]]);
    assert_eq!(
        bits(&solve_upper(&t, vecn([5.0, 4.0])).0),
        bits(&[1.0, 1.0])
    );
}

#[test]
fn a_zero_pivot_is_not_positive_definite() {
    // The boundary of the mask: `0 < d` is false at `d = +-0`, true at the smallest subnormal.
    for z in [0.0_f64, -0.0] {
        let (l, pd) = chol(&mat::<f64, 1>([[z]]));
        assert!(!pd && l.get(0, 0).to_bits() == 1.0_f64.to_bits());
    }
    let (l, pd) = chol(&mat::<f64, 1>([[f64::from_bits(1)]]));
    assert!(pd && l.get(0, 0) > 0.0 && l.get(0, 0).is_finite());
    // Exactly singular: `[[1, 1], [1, 1]]` has `d_2 = 1 - 1 = 0`; rank one `v v^T` for
    // `v = (1, 2, 3)` has `d_2 = 4 - 2^2 = 0` and `d_3 = 9 - 3^2 - 0 = 0`.
    let (l, pd) = chol(&mat::<f64, 2>([[1.0, 1.0], [1.0, 1.0]]));
    assert!(!pd && finite(&flat(&l)));
    let v = [[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [3.0, 6.0, 9.0]];
    let (l, pd) = chol(&mat::<f32, 3>(v));
    assert!(!pd && finite(&flat(&l)));
    let (l, pd) = chol(&mat::<D, 3>(v));
    assert!(!pd && finite(&flat(&l)));
    // A failed pivot gives `L_jj = 1`.
    assert_eq!(l.get(1, 1).v.to_bits(), 1.0_f64.to_bits());
    let (l, pd) = chol(&Mat3::<f64>::identity().scale(0.0));
    assert!(!pd && finite(&flat(&l)));
}

/// `v v^T + eps I` for `v = (1, 2, 3)` has pivots `1 + eps`, about `5 eps`, `~ eps`: set as soon as
/// `eps` clears the rounding of a pivot of size 14, `14 u`. Far above that is a set mask at both
/// widths, and the factor reads back.
#[test]
fn nearly_singular_is_positive_definite_above_the_rounding_of_a_pivot() {
    fn go<S: Lane>(eps: f64) -> Result<(), TestCaseError> {
        let v = [1.0, 2.0, 3.0];
        let a = mat::<S, 3>(array::from_fn(|r| {
            array::from_fn(|c| v[r] * v[c] + if r == c { eps } else { 0.0 })
        }));
        let (l, pd) = chol(&a);
        prop_assert!(pd.all(), "eps = {eps}");
        reconstructs(&a, &l)
    }
    for k in 0..=30 {
        let eps = libm::ldexp(1.0, -k);
        assert!(go::<f64>(eps).is_ok());
        if eps > 2048.0 * unit::<f32>() {
            assert!(go::<f32>(eps).is_ok());
            assert!(go::<D>(eps).is_ok());
        }
    }
}

fn clear<S: Lane>(c: [[f64; 2]; 2]) {
    let (l, pd) = chol(&mat::<S, 2>(c));
    assert!(!pd.any(), "{c:?}");
    assert!(finite(&flat(&l)), "{c:?}");
}

#[test]
fn indefinite_inputs_clear_the_mask() {
    let indefinite = [
        [[1.0, 2.0], [2.0, 1.0]],
        [[-1.0, 0.0], [0.0, 1.0]],
        [[1.0, 0.0], [0.0, -1.0]],
        [[0.0, 1.0], [1.0, 0.0]],
    ];
    for c in indefinite {
        clear::<f64>(c);
        clear::<f32>(c);
        clear::<D>(c);
    }
}

/// A NaN or `inf` in the lower triangle clears the mask and leaves a finite factor: a bad pivot
/// gives `L_jj = 1`, a bad entry below it is stored as `+0`, so every case here factors to `I`.
#[test]
fn non_finite_input_clears_the_mask_and_leaves_a_finite_factor() {
    fn go<S: Lane>(c: [[f64; 2]; 2]) {
        let (l, pd) = chol(&mat::<S, 2>(c));
        assert!(!pd.any(), "{c:?}");
        assert_eq!(
            bits(&flat(&l)),
            bits(&flat(&Matrix::<S, 2, 2>::identity())),
            "{c:?}"
        );
    }
    let (n, i) = (f64::NAN, f64::INFINITY);
    for c in [
        [[n, 0.0], [0.0, 1.0]],
        [[1.0, 0.0], [0.0, n]],
        [[1.0, 0.0], [n, 1.0]],
        [[i, 0.0], [0.0, 1.0]],
        [[1.0, 0.0], [i, 1.0]],
        [[1.0, 0.0], [-i, 1.0]],
    ] {
        go::<f64>(c);
        go::<f32>(c);
        go::<D>(c);
    }
    // The unread upper triangle is not consulted at all.
    let (l, pd) = chol(&mat::<f64, 2>([[1.0, n], [0.0, 1.0]]));
    assert!(pd && finite(&flat(&l)));
}

/// `M + M^T` from a pool with entries `[3.7e-4, 19.9, 0.99, 16.6, 0.80, 0.51]` at flat positions
/// `0, 2, 3, 6, 24, 30`: column 1's pivot fails at about `-4e5`, and letting `L_21` feed the
/// following columns squares the growth per column, so `f32` overflowed at `N = 6`.
#[test]
fn a_failed_pivot_does_not_feed_later_columns() {
    let mut pl = [[0.0; 3]; 64];
    let at = [
        (0, 3.7258983733806113e-4),
        (2, 19.867457469394203),
        (3, 0.9861989917093157),
    ];
    for (i, v) in at.into_iter().chain([
        (6, 16.581433788577492),
        (24, 0.7998084108647727),
        (30, 0.5089659033780938),
    ]) {
        pl[i][0] = v;
    }
    assert!(symmetric::<f32, 6>(&pl).is_ok());
    assert!(symmetric::<f64, 6>(&pl).is_ok());
}

/// A failed pivot zeroes its column, so garbage does not compound (without it, `M + M^T` of
/// entries below 20 overflows `f32` at `N = 6`). A quotient that overflows under a positive pivot
/// is a failure too: the entry is stored as `+0` and the mask clears, whatever the pivot's size.
/// Each input below has a finite lower triangle and used to give `inf` (a tiny pivot under huge
/// entries; a product of two finite entries that overflows in the sum).
#[test]
fn an_overflowing_entry_clears_the_mask_and_is_stored_as_zero() {
    let p = |e| libm::ldexp(1.0, e);
    // `d_1 = 2^-1074` makes `L_21 = 2^300`, `L_31 = 2^800`, and `L_31 L_21 = 2^1100` is `inf`.
    let big = mat::<f64, 3>([
        [p(-1074), p(-237), p(263)],
        [p(-237), p(602), 0.0],
        [p(263), 0.0, p(1000)],
    ]);
    // A pivot of `1e-300` is not subnormal: `L_10 = 1e300 / 1e-150` is `inf`.
    let tiny = [[1e-300, 1e300], [1e300, 1.0]];
    // Every pivot is at least 1: `L_20 L_10 = 1e400` overflows the sum in `L_21`.
    let wide = [[1.0, 1e100, 1e300], [1e100, 1e201, 0.0], [1e300, 0.0, 1.0]];
    let f32_wide = [[1.0, 1e10, 1e30], [1e10, 1e21, 0.0], [1e30, 0.0, 1.0]];
    fn zero_and_clear<S: Lane, const N: usize>(a: Matrix<S, N, N>) -> Matrix<S, N, N> {
        let (l, pd) = chol(&a);
        assert!(!pd.any() && finite(&flat(&l)) && strictly_upper_is_plus_zero(&l));
        l
    }
    zero_and_clear(big);
    let l = zero_and_clear(mat::<f64, 2>(tiny));
    assert_eq!(l.get(1, 0).to_bits(), 0);
    assert_eq!(l.get(1, 1).to_bits(), 1.0_f64.to_bits());
    zero_and_clear(mat::<f64, 3>(wide));
    zero_and_clear(mat::<f32, 3>(f32_wide));
    zero_and_clear(mat::<D, 3>(wide));
    // With `d_1 = 2^-900` the entries are `L_21 = 2^213`, `L_31 = 2^713` and the product is
    // `2^926`: nothing overflows and the mask is the plain pivot test.
    let mut fine = mat::<f64, 3>([
        [p(-900), p(-237), p(263)],
        [p(-237), p(602), 0.0],
        [p(263), 0.0, p(1000)],
    ]);
    let (l, pd) = chol(&fine);
    assert!(!pd && finite(&flat(&l)));
    fine.set(2, 2, p(1000) * 4.0);
    assert!(finite(&flat(&chol(&fine).0)));
}

/// A failed pivot at `j` leaves `L_jj = 1` and the column below it exactly `+0`, at `N = 4` so
/// that both rows below it are checked. `d_1 = 1 - 2^2 < 0`.
#[test]
fn a_failed_pivot_gives_a_unit_diagonal_and_a_zero_column() {
    fn go<S: Lane>() {
        let (l, pd) = chol(&mat::<S, 4>([
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 1.0, 5.0, 6.0],
            [3.0, 5.0, 9.0, 7.0],
            [4.0, 6.0, 7.0, 20.0],
        ]));
        assert!(!pd.any());
        assert_eq!(l.get(1, 1).lane(0).to_bits(), 1.0_f64.to_bits());
        for r in 2..4 {
            assert_eq!(l.get(r, 1).lane(0).to_bits(), 0);
        }
        assert!(finite(&flat(&l)));
    }
    go::<f64>();
    go::<f32>();
    go::<D>();
}

/// Sums add left to right from the first term and the pivot is `a_jj` minus the finished sum, not a
/// running subtraction (D16). With squares `1e16, 1, 1` in one sum: left to right gives `1e16`
/// (each `+1` is a tie that rounds back to even), so `d_3 = (1e16 + 4) - 1e16 = 4` and `L_33 = 2`;
/// right to left, or a running subtraction, gives `d_3 = 2`. The same terms in `L_43`'s inner
/// product pin the off-diagonal sum: `a_43 - 1e16 = 4`, `L_43 = 4 / 2`. (`L`, hand-built: `L_30 =
/// L_40 = 1e8`, `L_31 = L_32 = L_41 = L_42 = 1`, unit diagonal above.)
#[test]
fn the_factor_sums_left_to_right_and_subtracts_the_finished_sum() {
    let big = 1e8;
    let a = mat::<f64, 5>([
        [1.0, 0.0, 0.0, big, big],
        [0.0, 1.0, 0.0, 1.0, 1.0],
        [0.0, 0.0, 1.0, 1.0, 1.0],
        [big, 1.0, 1.0, 1e16 + 4.0, 1e16 + 4.0],
        [big, 1.0, 1.0, 1e16 + 4.0, 1e17],
    ]);
    let (l, _) = chol(&a);
    assert_eq!(l.get(3, 3).to_bits(), 2.0_f64.to_bits());
    assert_eq!(l.get(4, 3).to_bits(), 2.0_f64.to_bits());
}

/// The triangular solves add left to right too: with `l_30 x_0 = 1e16` and two unit terms, `b_3 =
/// 1e16 + 4` leaves `x_3 = 4`; right to left leaves `2`. The upper solve reads `u_0k x_k` in
/// increasing `k`, so the same terms in a row pin it.
#[test]
fn the_solves_sum_left_to_right() {
    let big = 1e8;
    let l = mat::<f64, 4>([
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [big, 1.0, 1.0, 1.0],
    ]);
    let x = solve_lower(&l, vecn([big, 1.0, 1.0, 1e16 + 4.0]));
    assert_eq!(bits(&x.0), bits(&[big, 1.0, 1.0, 4.0]));
    let u = mat::<f64, 4>([
        [1.0, big, 1.0, 1.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]);
    let x = solve_upper(&u, vecn([1e16 + 4.0, big, 1.0, 1.0]));
    assert_eq!(bits(&x.0), bits(&[4.0, big, 1.0, 1.0]));
}

/// `l^T x = b` reads column `i` of `l` below its diagonal in increasing `k`, the transpose of
/// `the_solves_sum_left_to_right`'s upper case: `l_10 x_1 = 1e16` and two unit terms leave `x_0 = 4`.
#[test]
fn the_transposed_solve_sums_left_to_right() {
    let big = 1e8;
    let l = mat::<f64, 4>([
        [1.0, 0.0, 0.0, 0.0],
        [big, 1.0, 0.0, 0.0],
        [1.0, 0.0, 1.0, 0.0],
        [1.0, 0.0, 0.0, 1.0],
    ]);
    let b = vecn([1e16 + 4.0, big, 1.0, 1.0]);
    let x = solve_lower_transposed(&l, b);
    assert_eq!(bits(&x.0), bits(&[4.0, big, 1.0, 1.0]));
    assert_eq!(bits(&x.0), bits(&solve_upper(&l.transpose(), b).0));
}

/// The upper triangle of `l` is not read by `chol_solve`, and `N = 0` is the empty vector.
#[test]
fn chol_solve_reads_only_the_lower_triangle() {
    let l = mat::<f64, 3>([[2.0, 0.0, 0.0], [6.0, 1.0, 0.0], [-8.0, 5.0, 3.0]]);
    let mut junk = l;
    junk.set(0, 1, f64::NAN);
    junk.set(0, 2, f64::INFINITY);
    junk.set(1, 2, -0.0);
    let b = vecn([4.0, 12.0, -16.0]);
    assert_eq!(bits(&chol_solve(&l, b).0), bits(&chol_solve(&junk, b).0));
    let e = Matrix::<f64, 0, 0>::identity();
    assert!(chol_solve(&e, Vector([])).0.is_empty());
}

#[test]
fn edge_sizes_and_scales() {
    let (l, pd) = chol(&Matrix::<f64, 0, 0>::identity());
    assert!(pd && flat(&l).is_empty());
    let (l, pd) = chol(&mat::<f64, 1>([[9.0]]));
    assert!(pd && l.get(0, 0).to_bits() == 3.0_f64.to_bits());
    // Scale invariance: `c^2 A` factors to `c L` exactly for `c = 2^k`, away from the range's ends.
    let a = mat::<f64, 2>([[4.0, 2.0], [2.0, 5.0]]);
    for k in [-100, -1, 1, 100] {
        let c = libm::ldexp(1.0, k);
        let (l, pd) = chol(&a.scale(c * c));
        assert!(pd);
        assert_eq!(l.get(1, 1).to_bits(), (2.0 * c).to_bits());
        assert_eq!(l.get(1, 0).to_bits(), c.to_bits());
    }
}

/// Differentiating `L_11 = sqrt(a_11)`, `L_21 = a_21 / L_11`, `L_22 = sqrt(a_22 - L_21^2)` by hand
/// at `a = [[4, 2], [2, 5]]` with `a_11`, `a_21` the two variables: `L = [[2, 0], [1, 2]]`,
/// `dL_11 = (1/4, 0)`, `dL_21 = (-1/8, 1/2)`, `dL_22 = (1/16, -1/4)`, all exact in binary.
#[test]
fn the_dual_derivative_matches_a_hand_derivation() {
    let a = Matrix::from_rows([
        Vector([D::variable(4.0, 0), D::variable(2.0, 1)]),
        Vector([D::variable(2.0, 1), D::constant(5.0)]),
    ]);
    let (l, pd) = chol(&a);
    assert!(pd);
    let want = [
        (0, 0, [0.25, 0.0]),
        (1, 0, [-0.125, 0.5]),
        (1, 1, [0.0625, -0.25]),
    ];
    for (r, c, d) in want {
        assert_eq!(bits(&l.get(r, c).d), bits(&d), "d L_{r}{c}");
    }
}

/// A zero diagonal divides by zero and returns `inf`/NaN; it never panics in release. (Debug
/// builds assert the domain, below.)
#[cfg(not(debug_assertions))]
#[test]
fn a_zero_diagonal_does_not_panic_in_release() {
    let z = mat::<f64, 2>([[0.0, 1.0], [1.0, 1.0]]);
    let (lo, up) = (
        solve_lower(&z, vecn([1.0, 1.0])),
        solve_upper(&z, vecn([1.0, 1.0])),
    );
    core::hint::black_box((lo, up));
}

/// Zero and NaN diagonals give the reference's bits in release, where nothing is asserted.
#[cfg(not(debug_assertions))]
#[test]
fn chol_solve_matches_the_reference_out_of_domain_in_release() {
    let n = f64::NAN;
    for d in [0.0, -0.0, n, -n] {
        let l = mat::<f64, 3>([[2.0, 0.0, 0.0], [6.0, d, 0.0], [-8.0, 5.0, 3.0]]);
        let b = vecn([4.0, 12.0, -16.0]);
        let want = solve_upper(&l.transpose(), solve_lower(&l, b));
        assert_eq!(all_bits(&chol_solve(&l, b).0), all_bits(&want.0));
    }
}

#[cfg(debug_assertions)]
mod out_of_domain {
    use super::*;

    #[test]
    #[should_panic(expected = "solve_lower: zero diagonal")]
    fn solve_lower_asserts_a_nonzero_diagonal() {
        let z = mat::<f64, 2>([[1.0, 0.0], [1.0, 0.0]]);
        core::hint::black_box(solve_lower(&z, vecn([1.0, 1.0])));
    }

    #[test]
    #[should_panic(expected = "solve_upper: zero diagonal")]
    fn solve_upper_asserts_a_nonzero_diagonal() {
        let z = mat::<f64, 2>([[0.0, 1.0], [0.0, 1.0]]);
        core::hint::black_box(solve_upper(&z, vecn([1.0, 1.0])));
    }

    /// The forward solve runs first and sees the same diagonal, so it is the one that fires.
    #[test]
    #[should_panic(expected = "solve_lower: zero diagonal")]
    fn chol_solve_asserts_a_nonzero_diagonal() {
        let z = mat::<f64, 2>([[1.0, 0.0], [1.0, 0.0]]);
        core::hint::black_box(chol_solve(&z, vecn([1.0, 1.0])));
    }

    #[test]
    #[should_panic(expected = "solve_lower_transposed: zero diagonal")]
    fn the_transposed_solve_asserts_a_nonzero_diagonal() {
        let z = mat::<f64, 2>([[1.0, 0.0], [1.0, 0.0]]);
        core::hint::black_box(solve_lower_transposed(&z, vecn([1.0, 1.0])));
    }

    #[test]
    #[should_panic(expected = "solve_lower_transposed: zero diagonal")]
    fn the_transposed_solve_asserts_a_non_nan_diagonal() {
        let z = mat::<f64, 2>([[1.0, 0.0], [1.0, f64::NAN]]);
        core::hint::black_box(solve_lower_transposed(&z, vecn([1.0, 1.0])));
    }

    #[test]
    #[should_panic(expected = "solve_lower: zero diagonal")]
    fn a_nan_diagonal_is_out_of_domain_too() {
        let z = mat::<f64, 1>([[f64::NAN]]);
        core::hint::black_box(solve_lower(&z, vecn([1.0])));
    }
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn spd_matrices_reconstruct_from_their_factor(pl in pool()) {
        spd_all::<f64>(&pl)?;
        spd_all::<f32>(&pl)?;
        spd_all::<D>(&pl)?;
    }

    #[test]
    fn any_symmetric_matrix_gives_a_finite_backward_stable_factor(pl in pool()) {
        symmetric_all::<f64>(&pl)?;
        symmetric_all::<f32>(&pl)?;
        symmetric_all::<D>(&pl)?;
    }

    #[test]
    fn a_wide_dynamic_range_stays_finite(
        v in prop::collection::vec((-1.0f64..=1.0, -600i32..600), 36)
    ) {
        wide::<f64, 2>(&v, 1)?;
        wide::<f64, 3>(&v, 1)?;
        wide::<f64, 4>(&v, 1)?;
        wide::<f64, 6>(&v, 1)?;
        wide::<D, 4>(&v, 1)?;
        wide::<f32, 3>(&v, 8)?;
        wide::<f32, 6>(&v, 8)?;
    }

    #[test]
    fn the_dual_value_lane_factors_like_the_plain_matrix(pl in pool()) {
        dual_value_is_plain::<2>(&pl)?;
        dual_value_is_plain::<3>(&pl)?;
        dual_value_is_plain::<5>(&pl)?;
    }

    #[test]
    fn a_negative_diagonal_is_never_positive_definite(pl in pool()) {
        negative_diagonal::<f64, 4>(&pl)?;
        negative_diagonal::<f32, 4>(&pl)?;
        negative_diagonal::<D, 4>(&pl)?;
    }

    #[test]
    fn rank_deficient_to_well_posed_stays_finite_and_backward_stable(pl in pool()) {
        near_singular::<f64>(&pl)?;
        near_singular::<f32>(&pl)?;
        near_singular::<D>(&pl)?;
    }

    #[test]
    fn triangular_solves_leave_a_small_residual(pl in pool()) {
        triangular_all::<f64>(&pl)?;
        triangular_all::<f32>(&pl)?;
        triangular_all::<D>(&pl)?;
    }

    /// The reference twin (`NUMERICS.md` §14): `chol_solve` is bit-identical to the transposed-copy
    /// composition (NaN sign and payload aside), over every special value and all lanes.
    #[test]
    fn chol_solve_matches_reference(l in cells(36, 0..24), b in cells(6, 0..24)) {
        solve_bits_all::<f64>(&l, &b)?;
        solve_bits_all::<f32>(&l, &b)?;
        solve_bits_all::<D>(&l, &b)?;
    }

    /// The same on plain values, where the rounding of every product and sum is exercised
    /// instead of the NaN and `inf` paths.
    #[test]
    fn chol_solve_matches_reference_on_plain_values(l in cells(36, 9..24), b in cells(6, 9..24)) {
        solve_bits_all::<f64>(&l, &b)?;
        solve_bits_all::<f32>(&l, &b)?;
        solve_bits_all::<D>(&l, &b)?;
    }

    #[test]
    fn cholesky_solve_leaves_a_small_residual(pl in pool()) {
        solve_spd_all::<f64>(&pl)?;
        solve_spd_all::<f32>(&pl)?;
        solve_spd_all::<D>(&pl)?;
    }
}
