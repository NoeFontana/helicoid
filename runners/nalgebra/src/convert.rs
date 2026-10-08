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
//!
//! The corpus's matrices are exactly symmetric, so the triangle read changes nothing; the corpus
//! test asserts it. Every record is answered in binary64, `@f32` strata included: the protocol
//! answers every stratum of a file, and the harness scores this runner at `f64` only.

use std::collections::BTreeMap;

use nalgebra::{Cholesky, Quaternion, SMatrix, SVector, Schur, SymmetricEigen, Unit};

/// Named arrays: a record's `in` object (its `shape` as two numbers), or the `out` of the answer.
pub(crate) type Fields = BTreeMap<String, Vec<f64>>;

/// One function id's evaluation.
type Evaluate = fn(&Fields) -> Result<Fields, String>;

/// The function ids this runner answers, each with its evaluation: the one list, so an id cannot
/// be advertised without an answer or answered without being advertised.
const TABLE: [(&str, Evaluate); 7] = [
    ("eig3", eig3),
    ("chol_n3", chol::<3>),
    ("chol_n6", chol::<6>),
    ("chol_solve_n3", chol_solve::<3>),
    ("chol_solve_n6", chol_solve::<6>),
    ("solve_cubic", solve_cubic),
    ("quat_renormalize", quat_renormalize),
];

/// The iterations `SymmetricEigen` and `Schur` get. Their `new` is `try_new(m, f64::EPSILON, 0)`,
/// which never gives up: the same answer wherever it converges within this bound, and a NaN
/// answer instead of a hang where it does not (`x³`, whose companion is a nilpotent Jordan block).
/// LAPACK's `dlahqr` allows `30 · max(10, n)`; the slowest convergent cubic found, `x³ ∓ 1`, takes
/// 39.
const MAX_ITERATIONS: usize = 300;

/// The function ids [`answer`] evaluates, in table order.
pub(crate) fn supported() -> impl Iterator<Item = &'static str> {
    TABLE.iter().map(|(id, _)| *id)
}

/// `fn_id`'s answer to `input`, or `None` if this runner does not answer `fn_id`.
pub(crate) fn answer(fn_id: &str, input: &Fields) -> Option<Result<Fields, String>> {
    TABLE
        .iter()
        .find(|(id, _)| *id == fn_id)
        .map(|(_, evaluate)| evaluate(input))
}

fn take<'a>(input: &'a Fields, key: &str, len: usize) -> Result<&'a [f64], String> {
    match input.get(key) {
        Some(v) if v.len() == len => Ok(v),
        Some(v) => Err(format!("`{key}` has {} values, expected {len}", v.len())),
        None => Err(format!("no input `{key}`")),
    }
}

/// The input `A`, which must carry `shape` `[N, N]`.
fn take_square<const N: usize>(input: &Fields) -> Result<&[f64], String> {
    let n = N as f64;
    match take(input, "shape", 2)? {
        [r, c] if r.to_bits() == n.to_bits() && c.to_bits() == n.to_bits() => {
            take(input, "A", N * N)
        }
        shape => Err(format!("`A` has shape {shape:?}, expected [{N}, {N}]")),
    }
}

/// A helicoid matrix, column-major, as nalgebra's: the same order.
pub(crate) fn helicoid_to_nalgebra_matrix<const R: usize, const C: usize>(
    m: &[f64],
) -> Result<SMatrix<f64, R, C>, String> {
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
pub(crate) fn nalgebra_to_helicoid_matrix<const R: usize, const C: usize>(
    m: &SMatrix<f64, R, C>,
) -> Vec<f64> {
    m.as_slice().to_vec()
}

/// A helicoid quaternion `[w, x, y, z]` as nalgebra's, whose `coords` are `[i, j, k, w]`.
pub(crate) fn helicoid_to_nalgebra_quat(q: &[f64]) -> Result<Quaternion<f64>, String> {
    match *q {
        [w, x, y, z] => Ok(Quaternion::new(w, x, y, z)),
        _ => Err(format!("a quaternion has 4 values, not {}", q.len())),
    }
}

/// A nalgebra quaternion as helicoid's `[w, x, y, z]`.
pub(crate) fn nalgebra_to_helicoid_quat(q: &Quaternion<f64>) -> Vec<f64> {
    vec![q.w, q.i, q.j, q.k]
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
fn eig3(input: &Fields) -> Result<Fields, String> {
    let a = helicoid_to_nalgebra_matrix::<3, 3>(take_square::<3>(input)?)?;
    let Some(eig) = SymmetricEigen::try_new(a, f64::EPSILON, MAX_ITERATIONS) else {
        return Ok(fields([
            ("V", vec![f64::NAN; 9]),
            ("lambda", vec![f64::NAN; 3]),
        ]));
    };
    let mut order = [0, 1, 2];
    order.sort_by(|&i, &j| eig.eigenvalues[i].total_cmp(&eig.eigenvalues[j]));
    let lambda = order.iter().map(|&i| eig.eigenvalues[i]).collect();
    let v = order
        .iter()
        .flat_map(|&i| {
            eig.eigenvectors
                .column(i)
                .iter()
                .copied()
                .collect::<Vec<_>>()
        })
        .collect();
    Ok(fields([("V", v), ("lambda", lambda)]))
}

/// `Cholesky::new(A)`: `valid = 1` and its `L` (zero above the diagonal), or `valid = 0` and zeros.
fn chol<const N: usize>(input: &Fields) -> Result<Fields, String> {
    let a = helicoid_to_nalgebra_matrix::<N, N>(take_square::<N>(input)?)?;
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
fn chol_solve<const N: usize>(input: &Fields) -> Result<Fields, String> {
    let a = helicoid_to_nalgebra_matrix::<N, N>(take_square::<N>(input)?)?;
    let b = SVector::<f64, N>::from_column_slice(take(input, "b", N)?);
    let x = Cholesky::new(a).map_or([f64::NAN; N].into(), |c| c.solve(&b));
    Ok(fields([("x", x.as_slice().to_vec())]))
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
fn solve_cubic(input: &Fields) -> Result<Fields, String> {
    let [a, b, c, d] = ["a", "b", "c", "d"].map(|k| take(input, k, 1).map(|v| v[0]));
    let (a, b, c, d) = (a?, b?, c?, d?);
    let none = || fields([("roots", vec![0.0; 3]), ("valid", vec![0.0; 3])]);
    let row = [-b / a, -c / a, -d / a];
    if a == 0.0 || ![a, b, c, d].iter().chain(&row).all(|x| x.is_finite()) {
        return Ok(none());
    }
    #[rustfmt::skip]
    let companion = SMatrix::<f64, 3, 3>::new(
        row[0], row[1], row[2],
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
    );
    let Some(schur) = Schur::try_new(companion, f64::EPSILON, MAX_ITERATIONS) else {
        return Ok(fields([
            ("roots", vec![f64::NAN; 3]),
            ("valid", vec![1.0; 3]),
        ]));
    };
    let (_, t) = schur.unpack();
    let (mut roots, mut valid) = (vec![0.0; 3], vec![0.0; 3]);
    let mut i = 0;
    while i < 3 {
        if i == 2 || t[(i + 1, i)] == 0.0 {
            (roots[i], valid[i]) = (t[(i, i)], 1.0);
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
fn quat_renormalize(input: &Fields) -> Result<Fields, String> {
    let mut q = Unit::new_unchecked(helicoid_to_nalgebra_quat(take(input, "q", 4)?)?);
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
        assert!(supported().all(|id| answer(id, &Fields::new()).is_some()));
        assert!(answer("so3_exp", &Fields::new()).is_none());
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
        let m = helicoid_to_nalgebra_matrix::<2, 3>(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0])?;
        assert_eq!(
            bits(&[m[(0, 0)], m[(0, 1)], m[(0, 2)]]),
            bits(&[1.0, 3.0, 5.0])
        );
        assert_eq!(m[(1, 0)].to_bits(), 2.0f64.to_bits());
        let back = nalgebra_to_helicoid_matrix(&m);
        assert_eq!(bits(&back), bits(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]));
        assert!(helicoid_to_nalgebra_matrix::<2, 2>(&[1.0; 3]).is_err());
        Ok(())
    }

    #[test]
    fn a_matrix_input_must_carry_its_shape() {
        let a = vec![0.0; 9];
        assert!(eig3(&input([("A", a.clone())])).is_err());
        assert!(eig3(&input([("A", a.clone()), ("shape", vec![3.0, 1.0])])).is_err());
        assert!(chol::<6>(&square::<3>(a)).is_err());
    }

    /// `diag(3, −1, 2)`: the eigenvalues ascending, `−1, 2, 3`, and the vectors the matching unit
    /// columns `e₁, e₂, e₀` up to sign.
    #[test]
    fn a_diagonal_eig3_is_sorted_ascending_with_its_unit_vectors() -> Result<(), String> {
        let out = eig3(&square::<3>(vec![
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
        let out = eig3(&square::<3>(vec![
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

    /// `[[4, 2, 2], [2, 5, 3], [2, 3, 6]] = L Lᵀ` with `L = [[2, 0, 0], [1, 2, 0], [1, 1, 2]]`:
    /// every step exact in binary64.
    #[test]
    fn an_spd_chol_is_valid_and_exact() -> Result<(), String> {
        let a = vec![4.0, 2.0, 2.0, 2.0, 5.0, 3.0, 2.0, 3.0, 6.0];
        let out = chol::<3>(&square::<3>(a.clone()));
        assert_eq!(bits(&get(&out, "valid")?), bits(&[1.0]));
        let l = [2.0, 1.0, 1.0, 0.0, 2.0, 1.0, 0.0, 0.0, 2.0];
        assert_eq!(bits(&get(&out, "L")?), bits(&l));
        // `A x = b` for `x = (1, −1, 2)`: `b = (6, 3, 11)`.
        let mut solve = square::<3>(a);
        solve.insert("b".into(), vec![6.0, 3.0, 11.0]);
        let x = get(&chol_solve::<3>(&solve), "x")?;
        for (got, want) in x.iter().zip([1.0, -1.0, 2.0]) {
            assert!((got - want).abs() < 1e-15, "{x:?}");
        }
        Ok(())
    }

    /// `[[1, 2, 0], [2, 1, 0], [0, 0, 1]]` has eigenvalues `−1, 1, 3`: the second pivot is `−3`.
    #[test]
    fn an_indefinite_chol_is_invalid_with_zeros_and_its_solve_is_nan() -> Result<(), String> {
        let a = vec![1.0, 2.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let out = chol::<3>(&square::<3>(a.clone()));
        assert_eq!(bits(&get(&out, "valid")?), bits(&[0.0]));
        assert_eq!(bits(&get(&out, "L")?), bits(&[0.0; 9]));
        let mut solve = square::<3>(a);
        solve.insert("b".into(), vec![1.0; 3]);
        assert!(get(&chol_solve::<3>(&solve), "x")?
            .iter()
            .all(|x| x.is_nan()));
        Ok(())
    }

    /// nalgebra reads the lower triangle: the upper one may hold anything.
    #[test]
    fn chol_reads_the_lower_triangle_only() -> Result<(), String> {
        let a = vec![4.0, 2.0, 2.0, 99.0, 5.0, 3.0, -7.0, f64::NAN, 6.0];
        let out = chol::<3>(&square::<3>(a));
        let l = [2.0, 1.0, 1.0, 0.0, 2.0, 1.0, 0.0, 0.0, 2.0];
        assert_eq!(bits(&get(&out, "L")?), bits(&l));
        Ok(())
    }

    fn cubic(a: f64, b: f64, c: f64, d: f64) -> Result<(Vec<f64>, Vec<f64>), String> {
        let out = solve_cubic(&input([
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
        let out = quat_renormalize(&input([("q", vec![w, 0.0, 0.0, 0.0])]));
        let want = 1.0 - 3.0 * 2f64.powi(-41);
        assert_eq!(bits(&get(&out, "q")?), bits(&[want, 0.0, 0.0, 0.0]));
        // Distinct components keep their places.
        let out = quat_renormalize(&input([("q", vec![0.5, -0.5, 0.25, 0.75])]));
        let f = 0.5 * (3.0 - 1.125);
        assert_eq!(
            bits(&get(&out, "q")?),
            bits(&[0.5 * f, -0.5 * f, 0.25 * f, 0.75 * f])
        );
        Ok(())
    }
}
