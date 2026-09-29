//! The corpus's arrays to sophus-rs's types and back, and the function ids this runner answers.
//!
//! `sophus_lie` 0.15.0 is rotation-first and `w`-first like `docs/decisions/0002`: a rotation is
//! `[re, im₀, im₁, im₂]`, an `Isometry3` is `[q; p]`, its tangent is `[ω; ν]`. Each conversion is a
//! layout change, not a permutation, and still a named function with a hand-computed test on values
//! whose components all differ: 0.10.0 to 0.14.0 are translation-first (`[υ; ω]`, `[p; q]`), so a
//! bump of the pin must fail here, not score (0014 (draft) question 28).
//! Inputs are handed over as the corpus has them, unnormalized. Where sophus-rs would `assert!`
//! the answer is an error: `|‖q‖ − 1| > 1e-6` (`q:nonunit` is at `2⁻⁴⁶`), or a rotation tangent
//! whose `‖φ‖²` is not finite, since `Exp` then builds a NaN quaternion and asserts on it.
//!
//! sophus-rs has the left Jacobian and its inverse; the right ones are answered at the negated
//! tangent, `J_r(τ) = J_l(−τ)` (`docs/NUMERICS.md` §1; question 29). Not answered: `coeff_*` (its
//! `a`, `b`, `c` are inline), `so3_act` (`transform`) and `sen3_ad_n1` (`adj`), both exported and
//! left to 0014 (draft) question 27, `so3_from_matrix`, and the `n2` and `n3` ids, since the crate
//! has no SE₂(3).

use std::collections::BTreeMap;

use sophus_autodiff::linalg::{MatF64, VecF64};
use sophus_lie::prelude::*;
use sophus_lie::{Isometry3F64, Rotation3F64};

/// Named arrays: a record's `in` object, or the `out` object of the answer.
pub(crate) type Fields = BTreeMap<String, Vec<f64>>;

/// One function id's evaluation.
type Evaluate = fn(&Fields) -> Result<Fields, String>;

/// Which Jacobian of the corpus an id asks for. sophus-rs has the left one and its inverse.
#[derive(Clone, Copy)]
enum Jacobian {
    Left,
    Right,
    LeftInverse,
    RightInverse,
}

impl Jacobian {
    /// The tangent to hand sophus-rs: `J_r(τ) = J_l(−τ)`, and likewise for the inverses.
    fn tangent<const N: usize>(self, tau: VecF64<N>) -> VecF64<N> {
        match self {
            Self::Left | Self::LeftInverse => tau,
            Self::Right | Self::RightInverse => -tau,
        }
    }

    fn is_inverse(self) -> bool {
        matches!(self, Self::LeftInverse | Self::RightInverse)
    }
}

/// The function ids this runner answers, each with its evaluation: the one list, so an id cannot
/// be advertised without an answer or answered without being advertised.
const TABLE: [(&str, Evaluate); 12] = [
    ("so3_exp", so3_exp),
    ("so3_log", so3_log),
    ("so3_jl", |i| so3_jacobian(i, Jacobian::Left)),
    ("so3_jr", |i| so3_jacobian(i, Jacobian::Right)),
    ("so3_jl_inv", |i| so3_jacobian(i, Jacobian::LeftInverse)),
    ("so3_jr_inv", |i| so3_jacobian(i, Jacobian::RightInverse)),
    ("sen3_exp_n1", sen3_exp_n1),
    ("sen3_log_n1", sen3_log_n1),
    ("sen3_jl_n1", |i| sen3_jacobian(i, Jacobian::Left)),
    ("sen3_jr_n1", |i| sen3_jacobian(i, Jacobian::Right)),
    ("sen3_jl_inv_n1", |i| {
        sen3_jacobian(i, Jacobian::LeftInverse)
    }),
    ("sen3_jr_inv_n1", |i| {
        sen3_jacobian(i, Jacobian::RightInverse)
    }),
];

/// The function ids [`answer`] evaluates, in table order.
pub(crate) fn supported() -> impl Iterator<Item = &'static str> {
    TABLE.iter().map(|(id, _)| *id)
}

fn take<'a>(input: &'a Fields, key: &str, len: usize) -> Result<&'a [f64], String> {
    match input.get(key) {
        Some(v) if v.len() == len => Ok(v),
        Some(v) => Err(format!("`{key}` has {} values, expected {len}", v.len())),
        None => Err(format!("no input `{key}`")),
    }
}

/// A helicoid quaternion `[w, x, y, z]` as sophus-rs's `[re, im₀, im₁, im₂]`: the same order.
pub(crate) fn helicoid_to_sophus_rs_quat(q: &[f64]) -> Result<VecF64<4>, String> {
    <[f64; 4]>::try_from(q)
        .map(VecF64::from)
        .map_err(|_| format!("a quaternion has 4 values, not {}", q.len()))
}

/// sophus-rs's `[re, im₀, im₁, im₂]` as a helicoid quaternion `[w, x, y, z]`.
pub(crate) fn sophus_rs_to_helicoid_quat(q: &VecF64<4>) -> Vec<f64> {
    q.as_slice().to_vec()
}

/// A helicoid 3-vector `[x, y, z]` as sophus-rs's: the same order.
pub(crate) fn helicoid_to_sophus_rs_vec3(v: &[f64]) -> Result<VecF64<3>, String> {
    <[f64; 3]>::try_from(v)
        .map(VecF64::from)
        .map_err(|_| format!("a vector has 3 values, not {}", v.len()))
}

/// sophus-rs's 3-vector as a helicoid `[x, y, z]`.
pub(crate) fn sophus_rs_to_helicoid_vec3(v: &VecF64<3>) -> Vec<f64> {
    v.as_slice().to_vec()
}

/// A helicoid SE(3) tangent `[φ; ρ]` as sophus-rs 0.15.0's `[ω; ν]`: the same order.
pub(crate) fn helicoid_to_sophus_rs_tangent(tau: &[f64]) -> Result<VecF64<6>, String> {
    <[f64; 6]>::try_from(tau)
        .map(VecF64::from)
        .map_err(|_| format!("a tangent has 6 values, not {}", tau.len()))
}

/// sophus-rs 0.15.0's `[ω; ν]` as a helicoid SE(3) tangent `[φ; ρ]`.
pub(crate) fn sophus_rs_to_helicoid_tangent(xi: &VecF64<6>) -> Vec<f64> {
    xi.as_slice().to_vec()
}

/// A helicoid `(q, x)` as an `Isometry3`'s parameters `[q; p]`.
pub(crate) fn helicoid_to_sophus_rs_iso3(q: &[f64], x: &[f64]) -> Result<VecF64<7>, String> {
    let (q, x) = (
        helicoid_to_sophus_rs_quat(q)?,
        helicoid_to_sophus_rs_vec3(x)?,
    );
    Ok(VecF64::from([q[0], q[1], q[2], q[3], x[0], x[1], x[2]]))
}

/// An `Isometry3`'s parameters `[q; p]` as a helicoid `(q, x)`.
pub(crate) fn sophus_rs_to_helicoid_iso3(params: &VecF64<7>) -> (Vec<f64>, Vec<f64>) {
    let (q, x) = params.as_slice().split_at(4);
    (q.to_vec(), x.to_vec())
}

/// A sophus-rs (nalgebra) matrix as the corpus's dense array: column-major, in the tangent's
/// order, which for `Isometry3`'s `[ω; ν]` is helicoid's rotation-first.
pub(crate) fn sophus_rs_to_helicoid_matrix<const R: usize, const C: usize>(
    m: &MatF64<R, C>,
) -> Vec<f64> {
    m.as_slice().to_vec()
}

fn fields<const N: usize>(entries: [(&str, Vec<f64>); N]) -> Fields {
    entries
        .into_iter()
        .map(|(k, v)| (k.to_string(), v))
        .collect()
}

/// `fn_id` evaluated at `input` by sophus-rs; `None` for an id this runner does not answer.
pub(crate) fn answer(fn_id: &str, input: &Fields) -> Option<Result<Fields, String>> {
    let (_, evaluate) = TABLE.iter().find(|(id, _)| *id == fn_id)?;
    Some(evaluate(input))
}

/// The rotation of `q`, or an error where sophus-rs's `from_params` would `assert!`.
fn rotation(q: &[f64]) -> Result<Rotation3F64, String> {
    let params = helicoid_to_sophus_rs_quat(q)?;
    match Rotation3F64::are_params_valid(params) {
        true => Ok(Rotation3F64::from_params(params)),
        false => Err(format!("sophus_lie refuses q = {q:?}: |‖q‖ − 1| > 1e-6")),
    }
}

/// The isometry of `(q, x)`, or an error where sophus-rs's `from_params` would `assert!`.
fn isometry(q: &[f64], x: &[f64]) -> Result<Isometry3F64, String> {
    let params = helicoid_to_sophus_rs_iso3(q, x)?;
    match Isometry3F64::are_params_valid(params) {
        true => Ok(Isometry3F64::from_params(params)),
        false => Err(format!("sophus_lie refuses q = {q:?}: |‖q‖ − 1| > 1e-6")),
    }
}

/// `phi`, or an error where sophus-rs's `exp` would `assert!` on its own NaN output: `‖φ‖²`
/// overflows from `|φ| ≈ 1.34e154`, and is NaN for a NaN component.
fn exp_domain(phi: VecF64<3>) -> Result<VecF64<3>, String> {
    match phi.squared_norm().is_finite() {
        true => Ok(phi),
        false => Err(format!(
            "sophus_lie cannot Exp φ = {phi:?}: ‖φ‖² is not finite"
        )),
    }
}

fn so3_exp(input: &Fields) -> Result<Fields, String> {
    let phi = exp_domain(helicoid_to_sophus_rs_vec3(take(input, "phi", 3)?)?)?;
    let q = sophus_rs_to_helicoid_quat(Rotation3F64::exp(phi).params());
    Ok(fields([("q", q)]))
}

fn so3_log(input: &Fields) -> Result<Fields, String> {
    let phi = rotation(take(input, "q", 4)?)?.log();
    Ok(fields([("phi", sophus_rs_to_helicoid_vec3(&phi))]))
}

fn sen3_exp_n1(input: &Fields) -> Result<Fields, String> {
    let tau = helicoid_to_sophus_rs_tangent(take(input, "tau", 6)?)?;
    exp_domain(tau.fixed_rows::<3>(0).into_owned())?;
    let (q, x) = sophus_rs_to_helicoid_iso3(Isometry3F64::exp(tau).params());
    Ok(fields([("q", q), ("x", x)]))
}

fn sen3_log_n1(input: &Fields) -> Result<Fields, String> {
    let tau = isometry(take(input, "q", 4)?, take(input, "x", 3)?)?.log();
    Ok(fields([("tau", sophus_rs_to_helicoid_tangent(&tau))]))
}

fn so3_jacobian(input: &Fields, which: Jacobian) -> Result<Fields, String> {
    let phi = which.tangent(helicoid_to_sophus_rs_vec3(take(input, "phi", 3)?)?);
    let j = match which.is_inverse() {
        true => Rotation3F64::inv_left_jacobian(phi),
        false => Rotation3F64::left_jacobian(phi),
    };
    Ok(fields([("J", sophus_rs_to_helicoid_matrix(&j))]))
}

fn sen3_jacobian(input: &Fields, which: Jacobian) -> Result<Fields, String> {
    let tau = which.tangent(helicoid_to_sophus_rs_tangent(take(input, "tau", 6)?)?);
    let j = match which.is_inverse() {
        true => Isometry3F64::inv_left_jacobian(tau),
        false => Isometry3F64::left_jacobian(tau),
    };
    Ok(fields([("J", sophus_rs_to_helicoid_matrix(&j))]))
}

#[cfg(test)]
mod tests {
    use core::f64::consts::{FRAC_1_SQRT_2, FRAC_2_PI, FRAC_PI_2, FRAC_PI_4};

    use super::*;

    /// Far above the rounding of the inputs and of sophus-rs at these angles, far below any
    /// convention slip.
    const TOL: f64 = 1e-14;

    /// `φ = (π/2)·(2, 3, 6)/7`: a quarter turn about an axis with rational, distinct components.
    fn skew_phi() -> [f64; 3] {
        [
            FRAC_PI_2 * 2.0 / 7.0,
            FRAC_PI_2 * 3.0 / 7.0,
            FRAC_PI_2 * 6.0 / 7.0,
        ]
    }

    /// Its unit quaternion, `(cos π/4, sin π/4 · (2, 3, 6)/7)`.
    fn skew_q() -> [f64; 4] {
        let s = FRAC_1_SQRT_2;
        [s, s * 2.0 / 7.0, s * 3.0 / 7.0, s * 6.0 / 7.0]
    }

    /// The `SKEW_*` constants are printed by `skew_constants.py` (`uv run`) from definitions at 60
    /// digits, no closed form of what the runner asks of sophus-rs. Here `x` of
    /// `Exp([skew_phi(); (1, 2, 3)])`, from `mpmath.expm` of the 4×4 matrix.
    const SKEW_X: [f64; 3] = [0.7494107237383275, 1.85168153974187, 3.157688988882956];

    /// `J_l(skew_phi())` and its inverse, row by row: `Σ Wⁿ/(n+1)!` and `mp.inverse`.
    const SKEW_JL: [[f64; 3]; 3] = [
        [0.6662834644192074, -0.501178552523345, 0.36182812145527005],
        [0.590169628678223, 0.7033630794837399, -0.04840474930127762],
        [-0.18384596914551396, 0.3153779777659117, 0.9035930008322155],
    ];
    const SKEW_JL_INV: [[f64; 3]; 3] = [
        [0.802916680671126, 0.6994762016797579, -0.28404366106358764],
        [-0.6469206498587249, 0.8248148272632231, 0.3032328029879634],
        [0.38915476470565374, -0.1455661475248642, 0.9430648188605475],
    ];
    /// The lower-left block of `J_l` and of `J_l⁻¹` at `τ = [skew_phi(); (1, 2, 3)]`: `ad_τ` from the
    /// commutator of the 4×4 hats, `Σ adⁿ/(n+1)!` and `mp.inverse`. The diagonal blocks are
    /// `SKEW_JL` and `SKEW_JL_INV`, the upper-right block is 0.
    const SKEW_Q: [[f64; 3]; 3] = [
        [-1.3905514865881037, -0.4568195866858147, 0.8780453915523152],
        [0.8672907514937674, -1.1479181324091472, 0.3950216908486655],
        [
            -0.18929451210682693,
            0.8363918035751929,
            -0.4722169806607858,
        ],
    ];
    const SKEW_Q_INV: [[f64; 3]; 3] = [
        [-0.9782953015843309, 1.642149639747513, -0.7547349411814258],
        [-1.357850360252487, -0.8175501960619139, 0.926448919242539],
        [1.245265058818574, -0.07355108075746103, -0.3242551447347997],
    ];

    fn near(got: &[f64], want: &[f64]) {
        assert_eq!(got.len(), want.len());
        for (g, w) in got.iter().zip(want) {
            assert!((g - w).abs() <= TOL, "{got:?} vs {want:?}");
        }
    }

    fn bits(v: &[f64]) -> Vec<u64> {
        v.iter().map(|x| x.to_bits()).collect()
    }

    fn input(entries: &[(&str, &[f64])]) -> Fields {
        entries
            .iter()
            .map(|(k, v)| (k.to_string(), v.to_vec()))
            .collect()
    }

    fn run(fn_id: &str, entries: &[(&str, &[f64])]) -> Result<Fields, String> {
        answer(fn_id, &input(entries)).ok_or_else(|| format!("{fn_id} is not answered"))?
    }

    /// The column-major dense array of rows `m`.
    fn dense3(m: [[f64; 3]; 3]) -> Vec<f64> {
        (0..3).flat_map(|c| m.map(|row| row[c])).collect()
    }

    fn transposed(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
        core::array::from_fn(|r| core::array::from_fn(|c| m[c][r]))
    }

    /// The column-major 6×6 of `[[D, 0], [Q, D]]`, blocks given by rows.
    fn dense6(d: [[f64; 3]; 3], q: [[f64; 3]; 3]) -> Vec<f64> {
        let entry = |r: usize, c: usize| match (r < 3, c < 3) {
            (true, true) | (false, false) => d[r % 3][c % 3],
            (false, true) => q[r - 3][c],
            (true, false) => 0.0,
        };
        (0..6)
            .flat_map(|c| (0..6).map(move |r| entry(r, c)))
            .collect()
    }

    /// `[ρ]×` for `ρ = (1, 2, 3)`, halved: the `Q` of `J_l` at `φ = 0`.
    const HALF_RHO_CROSS: [[f64; 3]; 3] = [[0.0, -1.5, 1.0], [1.5, 0.0, -0.5], [-1.0, 0.5, 0.0]];

    const EYE: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    #[test]
    fn a_quaternion_is_w_first_on_both_sides() -> Result<(), String> {
        // `q = (1, 2, 3, 4)/√30`; `R = (w² − |u|²)I + 2uuᵀ + 2w[u]×`, worked out by hand:
        // `30 R = [[-20, 4, 22], [20, -10, 20], [10, 28, 4]]`.
        let n = 30f64.sqrt();
        let q = [1.0 / n, 2.0 / n, 3.0 / n, 4.0 / n];
        let params = helicoid_to_sophus_rs_quat(&q)?;
        let r = Rotation3F64::from_params(params).matrix();
        let want = [[-20.0, 4.0, 22.0], [20.0, -10.0, 20.0], [10.0, 28.0, 4.0]];
        near(
            &sophus_rs_to_helicoid_matrix(&r),
            &dense3(want.map(|row| row.map(|e| e / 30.0))),
        );
        assert_eq!(bits(&sophus_rs_to_helicoid_quat(&params)), bits(&q));
        assert!(helicoid_to_sophus_rs_quat(&[1.0; 3]).is_err());
        Ok(())
    }

    #[test]
    fn a_vector_is_x_y_z_on_both_sides() -> Result<(), String> {
        let v = helicoid_to_sophus_rs_vec3(&[1.0, 2.0, 3.0])?;
        assert_eq!((v.x, v.y, v.z), (1.0, 2.0, 3.0));
        assert_eq!(
            sophus_rs_to_helicoid_vec3(&VecF64::<3>::new(4.0, 5.0, 6.0)),
            [4.0, 5.0, 6.0]
        );
        assert!(helicoid_to_sophus_rs_vec3(&[1.0; 2]).is_err());
        Ok(())
    }

    /// The hat of `[1, 2, 3, 4, 5, 6]` is `[[[ω]×, ν], [0, 0]]` with `ω = (1, 2, 3)` and
    /// `ν = (4, 5, 6)`: the tangent is rotation first on both sides. Translation-first, the block
    /// would be `[4, 5, 6]×` and `ν = (1, 2, 3)`.
    #[test]
    fn a_tangent_is_rotation_first_on_both_sides() -> Result<(), String> {
        let hat = Isometry3F64::hat(helicoid_to_sophus_rs_tangent(&[
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
        ])?);
        let rows = [
            [0.0, -3.0, 2.0, 4.0],
            [3.0, 0.0, -1.0, 5.0],
            [-2.0, 1.0, 0.0, 6.0],
            [0.0, 0.0, 0.0, 0.0],
        ];
        for (r, row) in rows.iter().enumerate() {
            assert_eq!(
                bits(&[hat[(r, 0)], hat[(r, 1)], hat[(r, 2)], hat[(r, 3)]]),
                bits(row)
            );
        }
        assert_eq!(
            sophus_rs_to_helicoid_tangent(&Isometry3F64::vee(hat)),
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        );
        assert!(helicoid_to_sophus_rs_tangent(&[0.0; 7]).is_err());
        Ok(())
    }

    /// `[R, x; 0, 1]` of `q = (1, 2, 3, 4)/√30` and `x = (5, 6, 7)`: the parameters are `[q; p]`.
    /// Translation-first, `x` would be read as `(w, i, j)`, `q` as `(k, p₀, p₁, p₂)`.
    #[test]
    fn an_isometry_is_q_then_x_on_both_sides() -> Result<(), String> {
        let n = 30f64.sqrt();
        let q = [1.0 / n, 2.0 / n, 3.0 / n, 4.0 / n];
        let params = helicoid_to_sophus_rs_iso3(&q, &[5.0, 6.0, 7.0])?;
        let m = Isometry3F64::from_params(params).matrix();
        let want = [
            [-20.0 / 30.0, 4.0 / 30.0, 22.0 / 30.0, 5.0],
            [20.0 / 30.0, -10.0 / 30.0, 20.0 / 30.0, 6.0],
            [10.0 / 30.0, 28.0 / 30.0, 4.0 / 30.0, 7.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        for (r, row) in want.iter().enumerate() {
            near(&[m[(r, 0)], m[(r, 1)], m[(r, 2)], m[(r, 3)]], row);
        }
        let (back_q, back_x) = sophus_rs_to_helicoid_iso3(&params);
        assert_eq!((bits(&back_q), back_x), (bits(&q), vec![5.0, 6.0, 7.0]));
        assert!(helicoid_to_sophus_rs_iso3(&[1.0; 4], &[1.0; 4]).is_err());
        assert!(helicoid_to_sophus_rs_iso3(&[1.0; 3], &[1.0; 3]).is_err());
        Ok(())
    }

    /// A dense matrix is column-major: the array of `[[1, 2, 3], [4, 5, 6], [7, 8, 9]]` is
    /// `[1, 4, 7, 2, 5, 8, 3, 6, 9]`, and of a 6×6, column j is the image of basis vector `e_j`.
    #[test]
    fn a_dense_matrix_is_column_major() {
        let m = MatF64::<3, 3>::new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
        assert_eq!(
            sophus_rs_to_helicoid_matrix(&m),
            [1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 9.0]
        );
        let mut six = MatF64::<6, 6>::zeros();
        six[(4, 1)] = 7.0;
        let dense = sophus_rs_to_helicoid_matrix(&six);
        // Entry (4, 1): column 1 starts at index 6.
        assert_eq!(
            (dense[6 + 4], dense.iter().filter(|x| **x != 0.0).count()),
            (7.0, 1)
        );
    }

    /// `Exp` of a quarter turn about `z`: `q = (cos π/4, 0, 0, sin π/4)`.
    #[test]
    fn so3_exp_of_a_quarter_turn_about_z() -> Result<(), String> {
        let out = run("so3_exp", &[("phi", &[0.0, 0.0, FRAC_PI_2])])?;
        near(&out["q"], &[FRAC_1_SQRT_2, 0.0, 0.0, FRAC_1_SQRT_2]);
        Ok(())
    }

    /// A rotation about `(2, 3, 6)/7` reads `x`, `y` and `z` apart; the quaternion is `w` first.
    #[test]
    fn so3_exp_about_a_skew_axis() -> Result<(), String> {
        near(&run("so3_exp", &[("phi", &skew_phi())])?["q"], &skew_q());
        Ok(())
    }

    /// `Log` reads `q` and `-q` as one rotation with `θ ∈ [0, π]`: sophus-rs flips on `w < 0`.
    #[test]
    fn so3_log_about_a_skew_axis_and_of_its_negative() -> Result<(), String> {
        let q = skew_q();
        for q in [q, q.map(|c| -c)] {
            near(&run("so3_log", &[("q", &q)])?["phi"], &skew_phi());
        }
        let s = FRAC_1_SQRT_2;
        near(
            &run("so3_log", &[("q", &[s, 0.0, 0.0, s])])?["phi"],
            &[0.0, 0.0, FRAC_PI_2],
        );
        Ok(())
    }

    /// `φ = 0`: `V = I`, so `x = ρ` and `q = 1`, whatever the layout of the six numbers.
    #[test]
    fn sen3_exp_and_log_of_a_pure_translation() -> Result<(), String> {
        let tau = [0.0, 0.0, 0.0, 1.0, 2.0, 3.0];
        let out = run("sen3_exp_n1", &[("tau", &tau)])?;
        near(&out["q"], &[1.0, 0.0, 0.0, 0.0]);
        near(&out["x"], &[1.0, 2.0, 3.0]);
        let back = run("sen3_log_n1", &[("q", &out["q"]), ("x", &out["x"])])?;
        near(&back["tau"], &tau);
        Ok(())
    }

    /// `τ = [φ; ρ] = [(π/2)(2, 3, 6)/7; (1, 2, 3)]`: the rotation is `skew_q()` and `x = V ρ` is
    /// `SKEW_X`. Read translation-first, the same six numbers would rotate by `(1, 2, 3)`.
    #[test]
    fn sen3_exp_about_a_skew_axis_with_a_skew_translation() -> Result<(), String> {
        let [a, b, c] = skew_phi();
        let out = run("sen3_exp_n1", &[("tau", &[a, b, c, 1.0, 2.0, 3.0])])?;
        near(&out["q"], &skew_q());
        near(&out["x"], &SKEW_X);
        Ok(())
    }

    #[test]
    fn sen3_log_inverts_it() -> Result<(), String> {
        let out = run("sen3_log_n1", &[("q", &skew_q()), ("x", &SKEW_X)])?;
        let [a, b, c] = skew_phi();
        near(&out["tau"], &[a, b, c, 1.0, 2.0, 3.0]);
        Ok(())
    }

    /// `x = V ρ` of the quarter turn about `z` and `ρ = (1, 0, 0)` is `(2/π, 2/π, 0)`.
    #[test]
    fn sen3_exp_of_a_quarter_turn_about_z_with_a_unit_translation() -> Result<(), String> {
        let tau = [0.0, 0.0, FRAC_PI_2, 1.0, 0.0, 0.0];
        let out = run("sen3_exp_n1", &[("tau", &tau)])?;
        near(&out["q"], &[FRAC_1_SQRT_2, 0.0, 0.0, FRAC_1_SQRT_2]);
        near(&out["x"], &[FRAC_2_PI, FRAC_2_PI, 0.0]);
        Ok(())
    }

    /// `W = (π/2)[e_z]×`, `W² = -(π/2)² diag(1, 1, 0)`: `J_l = I + (4/π²)W + ((4π−8)/π³)W²` is
    /// `[[2/π, −2/π, 0], [2/π, 2/π, 0], [0, 0, 1]]`, and `J_l⁻¹ = I − W/2 + (4/π² − 1/π)W²` is
    /// `[[π/4, π/4, 0], [−π/4, π/4, 0], [0, 0, 1]]`: `J_l J_l⁻¹ = I` by hand. `J_r(φ) = J_l(−φ)`
    /// is the transpose, `W` being skew and `W²` symmetric.
    #[test]
    fn so3_jacobians_of_a_quarter_turn_about_z() -> Result<(), String> {
        let (a, b) = (FRAC_2_PI, FRAC_PI_4);
        let jl = [[a, -a, 0.0], [a, a, 0.0], [0.0, 0.0, 1.0]];
        let jli = [[b, b, 0.0], [-b, b, 0.0], [0.0, 0.0, 1.0]];
        let phi = [0.0, 0.0, FRAC_PI_2];
        near(&run("so3_jl", &[("phi", &phi)])?["J"], &dense3(jl));
        near(&run("so3_jl_inv", &[("phi", &phi)])?["J"], &dense3(jli));
        near(
            &run("so3_jr", &[("phi", &phi)])?["J"],
            &dense3(transposed(jl)),
        );
        near(
            &run("so3_jr_inv", &[("phi", &phi)])?["J"],
            &dense3(transposed(jli)),
        );
        Ok(())
    }

    /// About a skew axis, and read row by row: an `x`/`y` swap of the input or a transposed output
    /// fails.
    #[test]
    fn so3_jacobians_about_a_skew_axis() -> Result<(), String> {
        let phi = skew_phi();
        near(&run("so3_jl", &[("phi", &phi)])?["J"], &dense3(SKEW_JL));
        near(
            &run("so3_jl_inv", &[("phi", &phi)])?["J"],
            &dense3(SKEW_JL_INV),
        );
        near(
            &run("so3_jr", &[("phi", &phi)])?["J"],
            &dense3(transposed(SKEW_JL)),
        );
        near(
            &run("so3_jr_inv", &[("phi", &phi)])?["J"],
            &dense3(transposed(SKEW_JL_INV)),
        );
        Ok(())
    }

    /// `φ = 0`: `ad_τ² = 0`, so `J_l = [[I, 0], [½[ρ]×, I]]` and `J_l⁻¹` has `−½[ρ]×`; `J_r(τ) =
    /// J_l(−τ)` swaps the signs. `ρ = (1, 2, 3)`: `[ρ]× = [[0, −3, 2], [3, 0, −1], [−2, 1, 0]]`.
    #[test]
    fn sen3_jacobians_of_a_pure_translation() -> Result<(), String> {
        let tau = [0.0, 0.0, 0.0, 1.0, 2.0, 3.0];
        let neg = HALF_RHO_CROSS.map(|row| row.map(|e| -e));
        let got = |id| Ok::<_, String>(run(id, &[("tau", &tau)])?["J"].clone());
        near(&got("sen3_jl_n1")?, &dense6(EYE, HALF_RHO_CROSS));
        near(&got("sen3_jl_inv_n1")?, &dense6(EYE, neg));
        near(&got("sen3_jr_n1")?, &dense6(EYE, neg));
        near(&got("sen3_jr_inv_n1")?, &dense6(EYE, HALF_RHO_CROSS));
        Ok(())
    }

    /// `τ = [skew_phi(); (1, 2, 3)]`: the diagonal blocks are `J_l(φ)`, the coupling block is `Q`.
    #[test]
    fn sen3_jacobians_about_a_skew_axis_with_a_skew_translation() -> Result<(), String> {
        let [a, b, c] = skew_phi();
        let tau = [a, b, c, 1.0, 2.0, 3.0];
        near(
            &run("sen3_jl_n1", &[("tau", &tau)])?["J"],
            &dense6(SKEW_JL, SKEW_Q),
        );
        near(
            &run("sen3_jl_inv_n1", &[("tau", &tau)])?["J"],
            &dense6(SKEW_JL_INV, SKEW_Q_INV),
        );
        Ok(())
    }

    /// The right Jacobians are the left ones at `−τ`, for the inverses too: a right one at `τ`
    /// is the left one at `−τ`, a left one is not the right one.
    #[test]
    fn a_right_jacobian_is_the_left_one_at_minus_tau() -> Result<(), String> {
        let [a, b, c] = skew_phi();
        let (tau, minus) = ([a, b, c, 1.0, 2.0, 3.0], [-a, -b, -c, -1.0, -2.0, -3.0]);
        for (right, left) in [
            ("sen3_jr_n1", "sen3_jl_n1"),
            ("sen3_jr_inv_n1", "sen3_jl_inv_n1"),
        ] {
            let (r, l) = (
                run(right, &[("tau", &tau)])?,
                run(left, &[("tau", &minus)])?,
            );
            assert_eq!(bits(&r["J"]), bits(&l["J"]));
            assert_ne!(bits(&r["J"]), bits(&run(left, &[("tau", &tau)])?["J"]));
        }
        Ok(())
    }

    /// The runner hands `q` over as the corpus has it: the answer is the library's on the raw
    /// input, not on its normalization. `q` is off the unit sphere by `1e-7`, far more than the
    /// corpus's `2⁻⁴⁶` and inside sophus-rs's `1e-6`, so that the two differ in the last bits
    /// (asserted).
    #[test]
    fn a_non_unit_quaternion_is_not_normalized() -> Result<(), String> {
        let q = skew_q().map(|c| c * (1.0 + 1e-7));
        let norm = q.iter().map(|c| c * c).sum::<f64>().sqrt();
        assert!((norm - 1.0 - 1e-7).abs() < 1e-9, "{norm}");
        let unit = q.map(|c| c / norm);
        let raw = run("so3_log", &[("q", &q)])?;
        assert_ne!(
            bits(&raw["phi"]),
            bits(&run("so3_log", &[("q", &unit)])?["phi"])
        );
        let want = Rotation3F64::from_params(VecF64::from(q)).log();
        assert_eq!(bits(&raw["phi"]), bits(want.as_slice()));
        let x = [1.0, 2.0, 3.0];
        let raw = run("sen3_log_n1", &[("q", &q), ("x", &x)])?;
        assert_ne!(
            bits(&raw["tau"]),
            bits(&run("sen3_log_n1", &[("q", &unit), ("x", &x)])?["tau"])
        );
        Ok(())
    }

    /// Where sophus-rs would `assert!`, the runner answers with an error naming the input: `q` off
    /// the unit sphere by more than `1e-6`, on either side, or not a number. Just inside: an answer.
    #[test]
    fn a_quaternion_sophus_rs_refuses_is_an_error_not_a_panic() {
        let x = [0.0; 3];
        for w in [
            2.0,
            0.5,
            1.0 + 2e-6,
            1.0 - 2e-6,
            0.0,
            f64::NAN,
            f64::INFINITY,
        ] {
            let q = [w, 0.0, 0.0, 0.0];
            assert!(run("so3_log", &[("q", &q)]).is_err(), "{w}");
            assert!(run("sen3_log_n1", &[("q", &q), ("x", &x)]).is_err(), "{w}");
        }
        for w in [1.0 + 9e-7, 1.0 - 9e-7] {
            let q = [w, 0.0, 0.0, 0.0];
            assert!(run("so3_log", &[("q", &q)]).is_ok(), "{w}");
            assert!(run("sen3_log_n1", &[("q", &q), ("x", &x)]).is_ok(), "{w}");
        }
    }

    /// `Exp` asserts on the NaN it makes of `‖φ‖² = ∞`, and only of that: `‖φ‖²` just below the
    /// overflow is answered. For SE(3) the rotation half decides, a huge `ρ` is an answer.
    #[test]
    fn a_tangent_whose_squared_norm_overflows_is_an_error_not_a_panic() {
        for big in [1.4e154, 1e300, f64::MAX, f64::INFINITY, f64::NAN] {
            assert!(
                run("so3_exp", &[("phi", &[big, 0.0, 0.0])]).is_err(),
                "{big}"
            );
            let tau = [0.0, big, 0.0, 1.0, 2.0, 3.0];
            assert!(run("sen3_exp_n1", &[("tau", &tau)]).is_err(), "{big}");
        }
        assert!(run("so3_exp", &[("phi", &[1e153, 1e153, 1e153])]).is_ok());
        let tau = [0.0, 0.0, 0.0, 1e300, 1e300, 1e300];
        assert!(run("sen3_exp_n1", &[("tau", &tau)]).is_ok());
    }

    /// No answered id panics on any of these values in any one place, or in all of them at once:
    /// an answer, finite or not, or an error. A panic fails the test.
    #[test]
    fn a_hostile_input_is_an_answer_or_an_error_never_a_panic() {
        let shape = |id: &str| -> &[(&'static str, usize)] {
            match id {
                "so3_log" => &[("q", 4)],
                "sen3_log_n1" => &[("q", 4), ("x", 3)],
                "so3_exp" | "so3_jl" | "so3_jr" | "so3_jl_inv" | "so3_jr_inv" => &[("phi", 3)],
                _ => &[("tau", 6)],
            }
        };
        let hostile = [
            0.0,
            -0.0,
            5e-324,
            1e-160,
            1.0,
            -1.0,
            core::f64::consts::PI,
            1.4e154,
            f64::MAX,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NAN,
        ];
        for id in supported() {
            let keys = shape(id);
            let len: usize = keys.iter().map(|(_, n)| n).sum();
            let split = |flat: &[f64]| -> Fields {
                let mut rest = flat;
                keys.iter()
                    .map(|(key, n)| {
                        let (head, tail) = rest.split_at(*n);
                        rest = tail;
                        (key.to_string(), head.to_vec())
                    })
                    .collect()
            };
            for v in hostile {
                let mut flats = vec![vec![v; len]];
                flats.extend((0..len).map(|at| {
                    let mut flat = vec![0.5; len];
                    flat[at] = v;
                    flat
                }));
                for flat in flats {
                    let _ = answer(id, &split(&flat));
                }
            }
        }
    }

    #[test]
    fn a_missing_or_short_input_is_an_error_and_an_unknown_id_is_none() {
        assert!(run("so3_exp", &[]).is_err());
        assert!(run("so3_exp", &[("phi", &[0.0; 4])]).is_err());
        assert!(run("sen3_jl_n1", &[("tau", &[0.0; 3])]).is_err());
        assert!(answer("coeff_a", &Fields::new()).is_none());
        assert!(answer("sen3_exp_n2", &Fields::new()).is_none());
    }

    /// The harness holds the same list (`RUNNERS` in `xtask`) and fails on a listed id that got no
    /// answer file, so dropping an id here is an error there, not a run with fewer rows.
    #[test]
    fn the_answered_ids_are_the_twelve_of_the_status_table() {
        let ids: Vec<&str> = supported().collect();
        assert_eq!(
            ids,
            [
                "so3_exp",
                "so3_log",
                "so3_jl",
                "so3_jr",
                "so3_jl_inv",
                "so3_jr_inv",
                "sen3_exp_n1",
                "sen3_log_n1",
                "sen3_jl_n1",
                "sen3_jr_n1",
                "sen3_jl_inv_n1",
                "sen3_jr_inv_n1"
            ]
        );
    }
}
