//! The corpus's arrays to nalgebra's types and back, and the function ids this runner answers
//! (`docs/decisions/0056` decision 4).
//!
//! A corpus matrix is column-major with a sibling `shape`, as nalgebra stores one, so a matrix
//! crosses unpermuted. A quaternion does not: helicoid's is `[w, x, y, z]` (`docs/decisions/0002`)
//! and nalgebra's `Quaternion::coords` is `[i, j, k, w]`, so it crosses through
//! `Quaternion::new(w, i, j, k)` and back through `.w`, `.i`, `.j`, `.k`.
//!
//! - `eig3`: `SymmetricEigen`, which reads the lower triangle only, eigenpairs sorted ascending.
//! - `chol_n3`, `chol_n6`: `Cholesky::new`, which reads the lower triangle only; `None` (a pivot
//!   that is zero, negative or NaN) is `valid = 0` with `L` zero.
//! - `chol_solve_n3`, `chol_solve_n6`: `Cholesky::new(A)` then `solve(b)`; `None` answers NaN, a
//!   non-finite row for the harness to record (`docs/PHASE1.md` §7), not an error.
//! - `solve_cubic`: the real Schur form of the monic companion matrix, `numpy.roots`' route without
//!   its balancing.
//! - `quat_renormalize`: `Unit::renormalize_fast`, the Newton step `q (3 − ‖q‖²)/2`.
//! - `mat2_inverse_adj` (`docs/decisions/0061`): `Matrix2::try_inverse`, what locus-tag calls, and
//!   `determinant`; `None` (a zero determinant) answers a NaN inverse, the singular answer.
//! - `svd3` (`docs/decisions/0101` decision 8): `SVD::try_new`, which nalgebra 0.35.0
//!   computes for a `Matrix3` by McAdams et al.'s route (`SymmetricEigen` of `AᵀA`, then a QR of
//!   `AV`), sorted descending and put in that record's signed convention.
//! - `nearest_rotation` (`docs/decisions/0102` decision 8, draft): `Rotation3::from_matrix_eps`,
//!   Müller et al.'s iteration, from the identity, `eps` the precision's unit roundoff, at most
//!   [`ROTATION_ITERATIONS`] steps.
//!
//! The matrices of `eig3` and `chol_*` are exactly symmetric, so the triangle read changes nothing;
//! the corpus test asserts it. In binary64 every record is answered, `@f32` strata included: the protocol
//! answers every stratum of a file. In binary32 (`docs/decisions/0077`) each id is the same
//! generic code at `f32`, so the iteration tolerance is `f32::EPSILON` where it was
//! `f64::EPSILON`; an input crosses by a lossless narrowing and an output by the exact widening.

use std::collections::BTreeMap;

use nalgebra::{
    Cholesky, Quaternion, RealField, Rotation3, SMatrix, SVector, Schur, SymmetricEigen, Unit, SVD,
};

/// Named arrays: a record's `in` object (its `shape` as two numbers), or the `out` of the answer.
pub(crate) type Fields = BTreeMap<String, Vec<f64>>;

/// The precision a record is answered in.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Precision {
    F64,
    F32,
}

/// The scalar an id is computed in. The corpus's values are binary64, so an `f32` input is the
/// binary32 it equals and an `f32` output the binary64 that equals it.
pub(crate) trait Real: RealField + Copy {
    /// `x` as this scalar, or `None` when that would round it.
    fn narrow(x: f64) -> Option<Self>;
    /// The binary64 equal to `self`.
    fn widen(self) -> f64;
}

impl Real for f64 {
    fn narrow(x: f64) -> Option<Self> {
        Some(x)
    }

    fn widen(self) -> f64 {
        self
    }
}

impl Real for f32 {
    fn narrow(x: f64) -> Option<Self> {
        // `as` rounds to nearest; the round trip says whether it had to.
        let y = x as f32;
        (f64::from(y).to_bits() == x.to_bits() || x.is_nan()).then_some(y)
    }

    fn widen(self) -> f64 {
        f64::from(self)
    }
}

/// One function id's evaluation.
type Evaluate = fn(&Fields) -> Result<Fields, String>;

/// The function ids this runner answers, each with its evaluation in binary64 and in binary32:
/// the one list, so an id cannot be advertised without an answer or answered without being
/// advertised.
const TABLE: [(&str, Evaluate, Evaluate); 10] = [
    ("eig3", eig3::<f64>, eig3::<f32>),
    ("chol_n3", chol::<f64, 3>, chol::<f32, 3>),
    ("chol_n6", chol::<f64, 6>, chol::<f32, 6>),
    ("chol_solve_n3", chol_solve::<f64, 3>, chol_solve::<f32, 3>),
    ("chol_solve_n6", chol_solve::<f64, 6>, chol_solve::<f32, 6>),
    ("solve_cubic", solve_cubic::<f64>, solve_cubic::<f32>),
    (
        "quat_renormalize",
        quat_renormalize::<f64>,
        quat_renormalize::<f32>,
    ),
    (
        "mat2_inverse_adj",
        mat2_inverse_adj::<f64>,
        mat2_inverse_adj::<f32>,
    ),
    ("svd3", svd3::<f64>, svd3::<f32>),
    (
        "nearest_rotation",
        nearest_rotation::<f64>,
        nearest_rotation::<f32>,
    ),
];

/// The iterations `SymmetricEigen`, `Schur` and `SVD` get. Their `new` is `try_new(m, T::EPSILON, 0)`,
/// which never gives up: the same answer wherever it converges within this bound, and a NaN
/// answer instead of a hang where it does not (`x³`, whose companion is a nilpotent Jordan block).
/// LAPACK's `dlahqr` allows `30 · max(10, n)`; the slowest convergent cubic found, `x³ ∓ 1`, takes
/// 39.
const MAX_ITERATIONS: usize = 300;

/// The steps `Rotation3::from_matrix_eps` gets (`docs/decisions/0102` decision 8, draft): at 100
/// it stops 0.80 short of the nearest rotation of `[[1,2,3],[4,5,6],[7,8,10]]`, at 1000 it is within
/// `5e-14`; its `0` would be unbounded.
const ROTATION_ITERATIONS: usize = 10_000;

/// The function ids [`answer`] evaluates, in table order.
pub(crate) fn supported() -> impl Iterator<Item = &'static str> {
    TABLE.iter().map(|(id, ..)| *id)
}

/// `fn_id`'s answer to `input` at `precision`, or `None` if this runner does not answer `fn_id`.
pub(crate) fn answer(
    fn_id: &str,
    input: &Fields,
    precision: Precision,
) -> Option<Result<Fields, String>> {
    let (_, wide, narrow) = TABLE.iter().find(|(id, ..)| *id == fn_id)?;
    Some(match precision {
        Precision::F64 => wide(input),
        Precision::F32 => narrow(input),
    })
}

fn take_f64<'a>(input: &'a Fields, key: &str, len: usize) -> Result<&'a [f64], String> {
    match input.get(key) {
        Some(v) if v.len() == len => Ok(v),
        Some(v) => Err(format!("`{key}` has {} values, expected {len}", v.len())),
        None => Err(format!("no input `{key}`")),
    }
}

/// The input `key`, `len` values, each exactly a `T`.
fn take<T: Real>(input: &Fields, key: &str, len: usize) -> Result<Vec<T>, String> {
    // Only a binary32 can refuse: no input is rounded on the way in.
    let narrow =
        |&x: &f64| T::narrow(x).ok_or_else(|| format!("`{key}` holds {x:e}, not a binary32"));
    take_f64(input, key, len)?.iter().map(narrow).collect()
}

fn take_one<T: Real>(input: &Fields, key: &str) -> Result<T, String> {
    Ok(take::<T>(input, key, 1)?[0])
}

/// The input `A`, which must carry `shape` `[N, N]`.
fn take_square<T: Real, const N: usize>(input: &Fields) -> Result<Vec<T>, String> {
    take_square_named::<T, N>(input, "A")
}

/// The input matrix `key`, which must carry `shape` `[N, N]`.
fn take_square_named<T: Real, const N: usize>(input: &Fields, key: &str) -> Result<Vec<T>, String> {
    let n = N as f64;
    match take_f64(input, "shape", 2)? {
        [r, c] if r.to_bits() == n.to_bits() && c.to_bits() == n.to_bits() => {
            take(input, key, N * N)
        }
        shape => Err(format!("`{key}` has shape {shape:?}, expected [{N}, {N}]")),
    }
}

/// `values` widened to binary64.
fn wide<T: Real>(values: &[T]) -> Vec<f64> {
    values.iter().map(|x| x.widen()).collect()
}

/// A helicoid matrix, column-major, as nalgebra's: the same order.
pub(crate) fn helicoid_to_nalgebra_matrix<T: Real, const R: usize, const C: usize>(
    m: &[T],
) -> Result<SMatrix<T, R, C>, String> {
    match m.len() == R * C {
        true => Ok(SMatrix::from_column_slice(m)),
        false => Err(format!(
            "a {R}x{C} matrix has {} values, not {}",
            R * C,
            m.len()
        )),
    }
}

/// A nalgebra matrix as helicoid's column-major array.
pub(crate) fn nalgebra_to_helicoid_matrix<T: Real, const R: usize, const C: usize>(
    m: &SMatrix<T, R, C>,
) -> Vec<f64> {
    wide(m.as_slice())
}

/// A helicoid quaternion `[w, x, y, z]` as nalgebra's, whose `coords` are `[i, j, k, w]`.
pub(crate) fn helicoid_to_nalgebra_quat<T: Real>(q: &[T]) -> Result<Quaternion<T>, String> {
    match *q {
        [w, x, y, z] => Ok(Quaternion::new(w, x, y, z)),
        _ => Err(format!("a quaternion has 4 values, not {}", q.len())),
    }
}

/// A nalgebra quaternion as helicoid's `[w, x, y, z]`.
pub(crate) fn nalgebra_to_helicoid_quat<T: Real>(q: &Quaternion<T>) -> Vec<f64> {
    wide(&[q.w, q.i, q.j, q.k])
}

fn fields<const N: usize>(entries: [(&str, Vec<f64>); N]) -> Fields {
    entries
        .into_iter()
        .map(|(k, v)| (k.to_string(), v))
        .collect()
}

fn flag(set: bool) -> f64 {
    if set {
        1.0
    } else {
        0.0
    }
}

/// `SymmetricEigen` of `A`, with the eigenpairs sorted ascending by eigenvalue (a stable sort, so
/// equal eigenvalues keep nalgebra's order) and each eigenvector column's sign as nalgebra gives
/// it. Not converged within [`MAX_ITERATIONS`]: NaN.
fn eig3<T: Real>(input: &Fields) -> Result<Fields, String> {
    let a = helicoid_to_nalgebra_matrix::<T, 3, 3>(&take_square::<T, 3>(input)?)?;
    let Some(eig) = SymmetricEigen::try_new(a, T::default_epsilon(), MAX_ITERATIONS) else {
        return Ok(fields([
            ("V", vec![f64::NAN; 9]),
            ("lambda", vec![f64::NAN; 3]),
        ]));
    };
    let values = wide(eig.eigenvalues.as_slice());
    let mut order = [0, 1, 2];
    order.sort_by(|&i, &j| values[i].total_cmp(&values[j]));
    let lambda = order.iter().map(|&i| values[i]).collect();
    let v = order
        .iter()
        .flat_map(|&i| {
            let column = eig.eigenvectors.column(i);
            column.iter().map(|x| x.widen()).collect::<Vec<_>>()
        })
        .collect();
    Ok(fields([("V", v), ("lambda", lambda)]))
}

/// `SVD::try_new(A, true, true, 5ε, MAX_ITERATIONS)`: `SVD::new`'s tolerance, with
/// [`MAX_ITERATIONS`] for its unbounded loop. Then `docs/decisions/0101` decision 1's convention:
/// the triples sorted descending by `σ` (a stable sort: McAdams' route orders by `‖AVeᵢ‖²` and
/// reports `|rᵢᵢ|`, which can invert a tie by an ulp); where `det V < 0`, `v₃` and `σ₃` negated; then where `det U < 0`, `u₃`
/// and `σ₃`. So `U, V ∈ SO(3)` and `σ₁ ≥ σ₂ ≥ |σ₃|`. Not converged: NaN.
fn svd3<T: Real>(input: &Fields) -> Result<Fields, String> {
    let a = helicoid_to_nalgebra_matrix::<T, 3, 3>(&take_square::<T, 3>(input)?)?;
    let tolerance = T::default_epsilon() * nalgebra::convert(5.0);
    let svd = SVD::try_new(a, true, true, tolerance, MAX_ITERATIONS);
    let Some((u, sigma, v_t)) = svd.and_then(|s| Some((s.u?, s.singular_values, s.v_t?))) else {
        return Ok(fields([
            ("U", vec![f64::NAN; 9]),
            ("V", vec![f64::NAN; 9]),
            ("sigma", vec![f64::NAN; 3]),
        ]));
    };
    let v = v_t.transpose();
    let mut order = [0, 1, 2];
    order.sort_by(|&i, &j| sigma[j].widen().total_cmp(&sigma[i].widen()));
    let mut u = SMatrix::<T, 3, 3>::from_columns(&order.map(|i| u.column(i).into_owned()));
    let mut v = SMatrix::<T, 3, 3>::from_columns(&order.map(|i| v.column(i).into_owned()));
    let mut sigma = order.map(|i| sigma[i]);
    for m in [&mut v, &mut u] {
        if m.determinant() < T::zero() {
            m.column_mut(2).neg_mut();
            sigma[2] = -sigma[2];
        }
    }
    Ok(fields([
        ("U", nalgebra_to_helicoid_matrix(&u)),
        ("V", nalgebra_to_helicoid_matrix(&v)),
        ("sigma", wide(&sigma)),
    ]))
}

/// `Rotation3::from_matrix_eps(M, u, ROTATION_ITERATIONS, I)`, `u = ε/2` the precision's unit
/// roundoff (`2⁻⁵³`, `2⁻²⁴`): the rotation nalgebra extracts from `M`, column-major.
fn nearest_rotation<T: Real>(input: &Fields) -> Result<Fields, String> {
    let m = helicoid_to_nalgebra_matrix::<T, 3, 3>(&take_square_named::<T, 3>(input, "M")?)?;
    let roundoff = T::default_epsilon() / nalgebra::convert(2.0);
    let r = Rotation3::from_matrix_eps(&m, roundoff, ROTATION_ITERATIONS, Rotation3::identity());
    Ok(fields([("R", nalgebra_to_helicoid_matrix(r.matrix()))]))
}

/// `A.try_inverse()` and `A.determinant()`; `None`: a NaN inverse.
fn mat2_inverse_adj<T: Real>(input: &Fields) -> Result<Fields, String> {
    let a = helicoid_to_nalgebra_matrix::<T, 2, 2>(&take_square::<T, 2>(input)?)?;
    let inv = a
        .try_inverse()
        .map_or(vec![f64::NAN; 4], |m| nalgebra_to_helicoid_matrix(&m));
    Ok(fields([
        ("det", vec![a.determinant().widen()]),
        ("inv", inv),
    ]))
}

/// `Cholesky::new(A)`: `valid = 1` and its `L` (zero above the diagonal), or `valid = 0` and zeros.
fn chol<T: Real, const N: usize>(input: &Fields) -> Result<Fields, String> {
    let a = helicoid_to_nalgebra_matrix::<T, N, N>(&take_square::<T, N>(input)?)?;
    let l = Cholesky::new(a).map(|c| c.l());
    Ok(fields([
        (
            "L",
            nalgebra_to_helicoid_matrix(&l.unwrap_or_else(SMatrix::zeros)),
        ),
        ("valid", vec![flag(l.is_some())]),
    ]))
}

/// `Cholesky::new(A)?.solve(b)`: forward then backward substitution. `None`: NaN.
fn chol_solve<T: Real, const N: usize>(input: &Fields) -> Result<Fields, String> {
    let a = helicoid_to_nalgebra_matrix::<T, N, N>(&take_square::<T, N>(input)?)?;
    let b = SVector::<T, N>::from_column_slice(&take::<T>(input, "b", N)?);
    let x = Cholesky::new(a).map_or(vec![f64::NAN; N], |c| wide(c.solve(&b).as_slice()));
    Ok(fields([("x", x)]))
}

/// The roots of `a x³ + b x² + c x + d` by `numpy.roots`' route, without LAPACK's balancing and
/// exceptional shifts: the eigenvalues of the monic
/// companion matrix
///
/// ```text
/// [ −b/a  −c/a  −d/a ]
/// [   1     0     0  ]
/// [   0     1     0  ]
/// ```
///
/// through its real Schur form `T` (`Schur::try_new`: Hessenberg reduction, then Francis's
/// implicit double-shift QR, with no balancing and no exceptional shift). A block of `T` ends where
/// its subdiagonal entry is exactly zero: a 1×1 block `T[i][i]` is a real root, a valid slot; a
/// 2×2 block is a complex pair, nalgebra having already split every 2×2 block whose eigenvalues
/// are real, so both its slots are invalid, value 0. The slots are in `T`'s order, not sorted.
/// `a = 0`, or a non-finite coefficient or companion entry: no valid slot. Not converged within
/// [`MAX_ITERATIONS`]: NaN in every slot, each marked valid, so the row is non-finite rather than
/// a scored miss.
fn solve_cubic<T: Real>(input: &Fields) -> Result<Fields, String> {
    let [a, b, c, d] = ["a", "b", "c", "d"].map(|k| take_one::<T>(input, k));
    let (a, b, c, d) = (a?, b?, c?, d?);
    let none = || fields([("roots", vec![0.0; 3]), ("valid", vec![0.0; 3])]);
    let row = [-b / a, -c / a, -d / a];
    if a == T::zero() || ![a, b, c, d].iter().chain(&row).all(|x| x.is_finite()) {
        return Ok(none());
    }
    let (one, zero) = (T::one(), T::zero());
    #[rustfmt::skip]
    let companion = SMatrix::<T, 3, 3>::new(
        row[0], row[1], row[2],
        one, zero, zero,
        zero, one, zero,
    );
    let Some(schur) = Schur::try_new(companion, T::default_epsilon(), MAX_ITERATIONS) else {
        return Ok(fields([
            ("roots", vec![f64::NAN; 3]),
            ("valid", vec![1.0; 3]),
        ]));
    };
    let (_, t) = schur.unpack();
    let (mut roots, mut valid) = (vec![0.0; 3], vec![0.0; 3]);
    let mut i = 0;
    while i < 3 {
        if i == 2 || t[(i + 1, i)] == zero {
            (roots[i], valid[i]) = (t[(i, i)].widen(), 1.0);
            i += 1;
        } else {
            i += 2;
        }
    }
    Ok(fields([("roots", roots), ("valid", valid)]))
}

/// `Unit::new_unchecked(q).renormalize_fast()`, which nalgebra 0.35.0 computes as
/// `s = (x² + z²) + (y² + w²)` (`dotx`'s four-lane pairing over `coords = [x, y, z, w]`), then
/// each component times `0.5 · (3 − s)`: the first-order Newton step of `NUMERICS.md` §3.6.
fn quat_renormalize<T: Real>(input: &Fields) -> Result<Fields, String> {
    let mut q = Unit::new_unchecked(helicoid_to_nalgebra_quat(&take::<T>(input, "q", 4)?)?);
    q.renormalize_fast();
    Ok(fields([("q", nalgebra_to_helicoid_quat(q.as_ref()))]))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn input<const N: usize>(entries: [(&str, Vec<f64>); N]) -> Fields {
        fields(entries)
    }

    fn square<const N: usize>(a: Vec<f64>) -> Fields {
        input([("A", a), ("shape", vec![N as f64, N as f64])])
    }

    fn get(out: &Result<Fields, String>, key: &str) -> Result<Vec<f64>, String> {
        let out = out.as_ref().map_err(Clone::clone)?;
        out.get(key).cloned().ok_or_else(|| format!("no `{key}`"))
    }

    fn bits(v: &[f64]) -> Vec<u64> {
        v.iter().map(|x| x.to_bits()).collect()
    }

    #[test]
    fn every_supported_id_has_an_answer_and_no_other() {
        assert_eq!(supported().count(), TABLE.len());
        assert!(supported().all(|id| answer(id, &Fields::new(), Precision::F64).is_some()));
        assert!(answer("so3_exp", &Fields::new(), Precision::F32).is_none());
    }

    /// Distinct components, so a permutation cannot pass.
    #[test]
    fn a_quaternion_crosses_w_first_into_ijkw_coords_and_back() -> Result<(), String> {
        let q = helicoid_to_nalgebra_quat(&[1.0, 2.0, 3.0, 4.0])?;
        assert_eq!(bits(q.coords.as_slice()), bits(&[2.0, 3.0, 4.0, 1.0]));
        assert_eq!(
            bits(&nalgebra_to_helicoid_quat(&q)),
            bits(&[1.0, 2.0, 3.0, 4.0])
        );
        assert!(helicoid_to_nalgebra_quat(&[1.0; 3]).is_err());
        Ok(())
    }

    /// `[1, 2, 3, 4, 5, 6]` column-major is the 2×3 matrix whose first row is `1, 3, 5`.
    #[test]
    fn a_matrix_crosses_column_major_both_ways() -> Result<(), String> {
        let m = helicoid_to_nalgebra_matrix::<f64, 2, 3>(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0])?;
        assert_eq!(
            bits(&[m[(0, 0)], m[(0, 1)], m[(0, 2)]]),
            bits(&[1.0, 3.0, 5.0])
        );
        assert_eq!(m[(1, 0)].to_bits(), 2.0f64.to_bits());
        let back = nalgebra_to_helicoid_matrix(&m);
        assert_eq!(bits(&back), bits(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]));
        assert!(helicoid_to_nalgebra_matrix::<f64, 2, 2>(&[1.0; 3]).is_err());
        Ok(())
    }

    #[test]
    fn a_matrix_input_must_carry_its_shape() {
        let a = vec![0.0; 9];
        assert!(eig3::<f64>(&input([("A", a.clone())])).is_err());
        assert!(eig3::<f64>(&input([("A", a.clone()), ("shape", vec![3.0, 1.0])])).is_err());
        assert!(svd3::<f64>(&input([("A", a.clone())])).is_err());
        assert!(nearest_rotation::<f64>(&input([("M", a.clone())])).is_err());
        assert!(nearest_rotation::<f64>(&square::<3>(a.clone())).is_err());
        assert!(chol::<f64, 6>(&square::<3>(a)).is_err());
    }

    /// `diag(3, −1, 2)`: the eigenvalues ascending, `−1, 2, 3`, and the vectors the matching unit
    /// columns `e₁, e₂, e₀` up to sign.
    #[test]
    fn a_diagonal_eig3_is_sorted_ascending_with_its_unit_vectors() -> Result<(), String> {
        let out = eig3::<f64>(&square::<3>(vec![
            3.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 2.0,
        ]));
        assert_eq!(bits(&get(&out, "lambda")?), bits(&[-1.0, 2.0, 3.0]));
        let v: Vec<f64> = get(&out, "V")?.iter().map(|x| x.abs()).collect();
        let want = [0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0];
        assert_eq!(bits(&v), bits(&want));
        Ok(())
    }

    /// `[[2, 1, 0], [1, 2, 0], [0, 0, 5]]`: eigenvalues `1, 3, 5`, the first two with vectors
    /// `(1, ∓1, 0)/√2`.
    #[test]
    fn a_coupled_eig3_has_its_closed_form_pairs() -> Result<(), String> {
        let out = eig3::<f64>(&square::<3>(vec![
            2.0, 1.0, 0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 5.0,
        ]));
        let lambda = get(&out, "lambda")?;
        for (got, want) in lambda.iter().zip([1.0, 3.0, 5.0]) {
            assert!((got - want).abs() < 1e-14, "{lambda:?}");
        }
        let v = get(&out, "V")?;
        let h = core::f64::consts::FRAC_1_SQRT_2;
        let col = |i: usize, want: [f64; 3]| {
            let c = &v[3 * i..3 * i + 3];
            let s = c[0].signum() * want[0].signum();
            let s = if want[0] == 0.0 { c[2].signum() } else { s };
            c.iter().zip(want).all(|(x, w)| (s * x - w).abs() < 1e-14)
        };
        assert!(
            col(0, [h, -h, 0.0]) && col(1, [h, h, 0.0]) && col(2, [0.0, 0.0, 1.0]),
            "{v:?}"
        );
        Ok(())
    }

    type M3 = SMatrix<f64, 3, 3>;

    /// The `svd3` answer as matrices: `U`, `σ`, `V`.
    fn svd3_of(a: [f64; 9]) -> Result<(M3, [f64; 3], M3), String> {
        let out = svd3::<f64>(&square::<3>(a.to_vec()));
        let u = helicoid_to_nalgebra_matrix::<f64, 3, 3>(&get(&out, "U")?)?;
        let v = helicoid_to_nalgebra_matrix::<f64, 3, 3>(&get(&out, "V")?)?;
        let sigma: [f64; 3] = get(&out, "sigma")?.try_into().map_err(|_| "sigma")?;
        Ok((u, sigma, v))
    }

    /// `U, V ∈ SO(3)` and `U diag(σ) Vᵀ = A`, each to `tol`.
    fn assert_signed_svd(a: [f64; 9], tol: f64) -> Result<[f64; 3], String> {
        let (u, sigma, v) = svd3_of(a)?;
        for (name, m) in [("U", u), ("V", v)] {
            assert!((m.determinant() - 1.0).abs() < tol, "det {name} {m}");
            let gram = m.transpose() * m - M3::identity();
            assert!(gram.norm() < tol, "{name}ᵀ{name} {gram}");
        }
        let s = M3::from_diagonal(&SVector::from(sigma));
        let back = u * s * v.transpose() - M3::from_column_slice(&a);
        assert!(back.norm() < tol, "U Σ Vᵀ − A = {back}");
        assert!(
            sigma[0] >= sigma[1] && sigma[1] >= sigma[2].abs(),
            "{sigma:?}"
        );
        Ok(sigma)
    }

    /// `diag(2, −3, 1)`, `det A < 0`: `σ = (3, 2, −1)`, the negative one the smallest in magnitude,
    /// whatever sign nalgebra's `U` and `V` came with.
    #[test]
    fn a_mixed_sign_diagonal_svd3_has_sigma3_negative_and_smallest() -> Result<(), String> {
        let sigma = assert_signed_svd([2.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 1.0], 1e-14)?;
        for (got, want) in sigma.iter().zip([3.0, 2.0, -1.0]) {
            assert!((got - want).abs() < 1e-14, "{sigma:?}");
        }
        Ok(())
    }

    /// `diag(1, 1, −5)`: `σ = (5, 1, −1)`, not `(1, 1, −5)`: `0015` PH.2's ordering, which
    /// `nearest_rotation = UVᵀ` needs.
    #[test]
    fn svd3_orders_by_magnitude_before_it_signs() -> Result<(), String> {
        let sigma = assert_signed_svd([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -5.0], 1e-14)?;
        for (got, want) in sigma.iter().zip([5.0, 1.0, -1.0]) {
            assert!((got - want).abs() < 1e-14, "{sigma:?}");
        }
        Ok(())
    }

    /// `R_z(π/2)`: `σ = (1, 1, 1)` and `UVᵀ = R`.
    #[test]
    fn a_rotation_svd3_has_unit_sigma_and_u_v_transpose_the_rotation() -> Result<(), String> {
        let r = [0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0];
        let sigma = assert_signed_svd(r, 1e-14)?;
        assert!(sigma.iter().all(|s| (s - 1.0).abs() < 1e-14), "{sigma:?}");
        let (u, _, v) = svd3_of(r)?;
        let diff = u * v.transpose() - M3::from_column_slice(&r);
        assert!(diff.norm() < 1e-14, "{diff}");
        Ok(())
    }

    /// General matrices of both signs of `det`: reconstruction, both determinants `+1`, the order;
    /// `sign σ₃ = sign det A`. In binary32 every output is a binary32.
    #[test]
    fn a_general_svd3_reconstructs_in_the_signed_convention() -> Result<(), String> {
        // `[[1, 2, 3], [4, 5, 6], [7, 8, 10]]`, `det = −3`, and its negation, `det = 3`.
        let negative = [1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 10.0];
        let positive = negative.map(|x| -x);
        for a in [positive, negative] {
            let det = M3::from_column_slice(&a).determinant();
            let sigma = assert_signed_svd(a, 1e-13)?;
            assert_eq!(sigma[2] < 0.0, det < 0.0, "{sigma:?} {det}");
        }
        let narrow = svd3::<f32>(&square::<3>(negative.to_vec()));
        for key in ["U", "V", "sigma"] {
            let values = get(&narrow, key)?;
            assert!(values
                .iter()
                .all(|&x| f64::from(x as f32).to_bits() == x.to_bits()));
        }
        let sigma = get(&narrow, "sigma")?;
        assert!(sigma[2] < 0.0 && sigma[1] >= sigma[2].abs(), "{sigma:?}");
        Ok(())
    }

    /// `nearest_rotation`'s `R` as a matrix.
    fn nearest_rotation_of<T: Real>(m: [f64; 9]) -> Result<M3, String> {
        let out = nearest_rotation::<T>(&input([("M", m.to_vec()), ("shape", vec![3.0, 3.0])]));
        helicoid_to_nalgebra_matrix::<f64, 3, 3>(&get(&out, "R")?)
    }

    /// A rotation about `(1, 2, 2)/3` by `0.7`: `R_z(π/2)` alone would be exact in every entry.
    fn turn() -> M3 {
        let axis = Unit::new_normalize(SVector::<f64, 3>::new(1.0, 2.0, 2.0));
        *Rotation3::from_axis_angle(&axis, 0.7).matrix()
    }

    /// An exact rotation comes back to rounding, `R_z(π/2)` bit for bit, column-major.
    #[test]
    fn nearest_rotation_returns_a_rotation_unchanged() -> Result<(), String> {
        let quarter = [0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0];
        let r = nearest_rotation_of::<f64>(quarter)?;
        assert!((r - M3::from_column_slice(&quarter)).norm() < 1e-15, "{r}");
        let t = turn();
        let r = nearest_rotation_of::<f64>(t.as_slice().try_into().map_err(|_| "9")?)?;
        assert!((r - t).norm() < 1e-14, "{r} {t}");
        Ok(())
    }

    /// `cR` for `c ∈ {10⁻³, 7, 10³}`: `R`, and in binary32 every output a binary32.
    #[test]
    fn nearest_rotation_of_a_scaled_rotation_is_the_rotation() -> Result<(), String> {
        let t = turn();
        for c in [1e-3, 7.0, 1e3] {
            let r = nearest_rotation_of::<f64>((t * c).as_slice().try_into().map_err(|_| "9")?)?;
            assert!((r - t).norm() < 1e-13, "c = {c}: {r}");
            assert!((r.determinant() - 1.0).abs() < 1e-13);
        }
        let narrow: [f64; 9] = t
            .map(|x| f64::from(x as f32) * 4.0)
            .as_slice()
            .try_into()
            .map_err(|_| "9")?;
        let r = nearest_rotation_of::<f32>(narrow)?;
        assert!(r
            .iter()
            .all(|&x| f64::from(x as f32).to_bits() == x.to_bits()));
        assert!((r - t).norm() < 1e-5, "{r}");
        Ok(())
    }

    /// `[[4, 2, 2], [2, 5, 3], [2, 3, 6]] = L Lᵀ` with `L = [[2, 0, 0], [1, 2, 0], [1, 1, 2]]`:
    /// every step exact in binary64.
    #[test]
    fn an_spd_chol_is_valid_and_exact() -> Result<(), String> {
        let a = vec![4.0, 2.0, 2.0, 2.0, 5.0, 3.0, 2.0, 3.0, 6.0];
        let out = chol::<f64, 3>(&square::<3>(a.clone()));
        assert_eq!(bits(&get(&out, "valid")?), bits(&[1.0]));
        let l = [2.0, 1.0, 1.0, 0.0, 2.0, 1.0, 0.0, 0.0, 2.0];
        assert_eq!(bits(&get(&out, "L")?), bits(&l));
        // `A x = b` for `x = (1, −1, 2)`: `b = (6, 3, 11)`.
        let mut solve = square::<3>(a);
        solve.insert("b".into(), vec![6.0, 3.0, 11.0]);
        let x = get(&chol_solve::<f64, 3>(&solve), "x")?;
        for (got, want) in x.iter().zip([1.0, -1.0, 2.0]) {
            assert!((got - want).abs() < 1e-15, "{x:?}");
        }
        Ok(())
    }

    /// `[[1, 2, 0], [2, 1, 0], [0, 0, 1]]` has eigenvalues `−1, 1, 3`: the second pivot is `−3`.
    #[test]
    fn an_indefinite_chol_is_invalid_with_zeros_and_its_solve_is_nan() -> Result<(), String> {
        let a = vec![1.0, 2.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let out = chol::<f64, 3>(&square::<3>(a.clone()));
        assert_eq!(bits(&get(&out, "valid")?), bits(&[0.0]));
        assert_eq!(bits(&get(&out, "L")?), bits(&[0.0; 9]));
        let mut solve = square::<3>(a);
        solve.insert("b".into(), vec![1.0; 3]);
        assert!(get(&chol_solve::<f64, 3>(&solve), "x")?
            .iter()
            .all(|x| x.is_nan()));
        Ok(())
    }

    /// nalgebra reads the lower triangle: the upper one may hold anything.
    #[test]
    fn chol_reads_the_lower_triangle_only() -> Result<(), String> {
        let a = vec![4.0, 2.0, 2.0, 99.0, 5.0, 3.0, -7.0, f64::NAN, 6.0];
        let out = chol::<f64, 3>(&square::<3>(a));
        let l = [2.0, 1.0, 1.0, 0.0, 2.0, 1.0, 0.0, 0.0, 2.0];
        assert_eq!(bits(&get(&out, "L")?), bits(&l));
        Ok(())
    }

    fn cubic(a: f64, b: f64, c: f64, d: f64) -> Result<(Vec<f64>, Vec<f64>), String> {
        let out = solve_cubic::<f64>(&input([
            ("a", vec![a]),
            ("b", vec![b]),
            ("c", vec![c]),
            ("d", vec![d]),
        ]));
        Ok((get(&out, "roots")?, get(&out, "valid")?))
    }

    /// `2 (x − 1)(x − 2)(x + 3) = 2x³ − 14x + 12`: three valid slots, `{−3, 1, 2}` in some order.
    #[test]
    fn a_cubic_with_three_real_roots_has_three_valid_slots() -> Result<(), String> {
        let (roots, valid) = cubic(2.0, 0.0, -14.0, 12.0)?;
        assert_eq!(bits(&valid), bits(&[1.0; 3]));
        let mut sorted = roots.clone();
        sorted.sort_by(f64::total_cmp);
        for (got, want) in sorted.iter().zip([-3.0, 1.0, 2.0]) {
            assert!((got - want).abs() < 1e-13, "{roots:?}");
        }
        Ok(())
    }

    /// `(x − 2)(x² + 1) = x³ − 2x² + x − 2`: one valid slot, `2`; the pair `±i` two invalid zeros.
    #[test]
    fn a_complex_pair_is_two_invalid_zero_slots() -> Result<(), String> {
        let (roots, valid) = cubic(1.0, -2.0, 1.0, -2.0)?;
        assert_eq!(valid.iter().sum::<f64>().to_bits(), 1.0f64.to_bits());
        for (r, v) in roots.iter().zip(&valid) {
            match v.to_bits() == 1.0f64.to_bits() {
                true => assert!((r - 2.0).abs() < 1e-14, "{roots:?}"),
                false => assert_eq!(r.to_bits(), 0.0f64.to_bits()),
            }
        }
        Ok(())
    }

    /// `x³ − 1`: the companion is a cyclic permutation, the case that stalls a double-shift QR
    /// without an exceptional shift; nalgebra converges, to the one real root `1`.
    #[test]
    fn the_cyclic_companion_converges_to_its_one_real_root() -> Result<(), String> {
        let (roots, valid) = cubic(1.0, 0.0, 0.0, -1.0)?;
        assert_eq!(valid.iter().sum::<f64>().to_bits(), 1.0f64.to_bits());
        let at = valid.iter().position(|v| v.to_bits() == 1.0f64.to_bits());
        assert!(
            at.is_some_and(|i| (roots[i] - 1.0).abs() < 1e-13),
            "{roots:?}"
        );
        Ok(())
    }

    /// `x³`: the companion is a nilpotent Jordan block, on which `Schur::new` never returns.
    #[test]
    fn a_cubic_whose_schur_form_does_not_converge_is_non_finite_not_a_hang() -> Result<(), String> {
        let (roots, valid) = cubic(1.0, 0.0, 0.0, 0.0)?;
        assert!(roots.iter().all(|r| r.is_nan()));
        assert_eq!(bits(&valid), bits(&[1.0; 3]));
        Ok(())
    }

    #[test]
    fn a_zero_or_non_finite_leading_coefficient_has_no_valid_slot() -> Result<(), String> {
        for a in [0.0, f64::NAN, f64::INFINITY] {
            let (roots, valid) = cubic(a, 1.0, 1.0, 1.0)?;
            assert_eq!(bits(&valid), bits(&[0.0; 3]));
            assert_eq!(bits(&roots), bits(&[0.0; 3]));
        }
        Ok(())
    }

    /// `q = (1 + 2⁻²⁰)(1, 0, 0, 0)`: `s = 1 + 2⁻¹⁹ + 2⁻⁴⁰` and `0.5 (3 − s) = 1 − 2⁻²⁰ − 2⁻⁴¹`
    /// exactly, and the product `1 − 3 · 2⁻⁴¹ − 2⁻⁶¹` rounds to `1 − 3 · 2⁻⁴¹`, all by hand: the
    /// Newton step, not the projection, which is `1`.
    #[test]
    fn renormalize_is_the_newton_step_in_w_first_order() -> Result<(), String> {
        let w = 1.0 + 2f64.powi(-20);
        let out = quat_renormalize::<f64>(&input([("q", vec![w, 0.0, 0.0, 0.0])]));
        let want = 1.0 - 3.0 * 2f64.powi(-41);
        assert_eq!(bits(&get(&out, "q")?), bits(&[want, 0.0, 0.0, 0.0]));
        // Distinct components keep their places.
        let out = quat_renormalize::<f64>(&input([("q", vec![0.5, -0.5, 0.25, 0.75])]));
        let f = 0.5 * (3.0 - 1.125);
        assert_eq!(
            bits(&get(&out, "q")?),
            bits(&[0.5 * f, -0.5 * f, 0.25 * f, 0.75 * f])
        );
        Ok(())
    }

    /// An input crosses into binary32 exactly or not at all, and an output back exactly.
    #[test]
    fn a_binary32_input_is_the_binary64_it_equals_and_no_other() {
        let one_ulp = 1.0 + 2f64.powi(-23);
        assert_eq!(
            f32::narrow(one_ulp).map(f32::to_bits),
            Some(1.0f32.next_up().to_bits())
        );
        assert_eq!(f32::narrow(1.0 + 2f64.powi(-30)), None);
        // A binary32 subnormal, and the binary64 nearest `10⁻⁴⁰`, which is none.
        let subnormal = f64::from(1e-40f32);
        let back = f32::narrow(subnormal).map(f32::widen);
        assert_eq!(back.map(f64::to_bits), Some(subnormal.to_bits()));
        assert_eq!(f32::narrow(1e-40), None);
        assert!(f32::narrow(f64::NAN).is_some_and(f32::is_nan));
        assert_eq!(f32::narrow(f64::MAX), None);
        let q = helicoid_to_nalgebra_quat::<f32>(&[1.0, 2.0, 3.0, 4.0]).map(|q| q.coords);
        assert_eq!(
            q.map(|c| c.as_slice().to_vec()),
            Ok(vec![2.0f32, 3.0, 4.0, 1.0])
        );
        let take = take::<f32>(&input([("q", vec![1.0 + 2f64.powi(-30)])]), "q", 1);
        assert!(take.is_err_and(|e| e.contains("not a binary32")));
    }

    /// `q = (1 + 2⁻¹⁰)(1, 0, 0, 0)` in binary32: `s = 1 + 2⁻⁹ + 2⁻²⁰`, `0.5 (3 − s) = 1 − 2⁻¹⁰ −
    /// 2⁻²¹` and the product `1 − 3 · 2⁻²¹ − 2⁻³¹` exactly, by hand, which binary32 rounds to
    /// `1 − 3 · 2⁻²¹` and binary64 keeps: the step is computed in binary32, not promoted.
    #[test]
    fn renormalize_in_binary32_rounds_in_binary32() -> Result<(), String> {
        let w = 1.0 + 2f64.powi(-10);
        let q = |out: &Result<Fields, String>| get(out, "q").map(|q| q[0]);
        let narrow = quat_renormalize::<f32>(&input([("q", vec![w, 0.0, 0.0, 0.0])]));
        let wide = quat_renormalize::<f64>(&input([("q", vec![w, 0.0, 0.0, 0.0])]));
        let rounded = 1.0 - 3.0 * 2f64.powi(-21);
        assert_eq!(q(&narrow)?.to_bits(), rounded.to_bits());
        assert_eq!(q(&wide)?.to_bits(), (rounded - 2f64.powi(-31)).to_bits());
        Ok(())
    }
}
