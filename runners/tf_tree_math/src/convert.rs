//! The corpus's arrays to `tf_tree_math`'s types and back, and the function ids this runner answers.
//!
//! `tf_tree_math`'s conventions are the ones `docs/decisions/0002` adopts: Hamilton `w` first,
//! twists `[ω, v]` (rotation first). So each conversion is a layout change, not a permutation; it
//! is still a named function with a hand-computed test, on values whose components all differ,
//! because that is what keeps a later change of either convention from passing silently
//! (`docs/PHASE1.md` §7).
//! Inputs are handed over as the corpus has them: no normalization, no clamping. `tf_tree_math`'s
//! own contract admits `|‖q‖² − 1| ≤ 1e-12`, which `q:nonunit` (`2⁻⁴⁵`) respects.
//!
//! Not answered, because `tf_tree_math` does not export them at the pin: `coeff_*` (`v_coeffs`
//! and `vinv_c3` are private) and every Jacobian `J` of `Exp` (there is no `J_r` or `J_l`).

use std::collections::BTreeMap;

use tf_tree_math::{exp_se3, exp_so3, log_se3, log_so3, Iso3, Quat, Vec3};

/// Named arrays: a record's `in` object, or the `out` object of the answer.
pub(crate) type Fields = BTreeMap<String, Vec<f64>>;

/// One function id's evaluation.
type Evaluate = fn(&Fields) -> Result<Fields, String>;

/// The function ids this runner answers, each with its evaluation: the one list, so an id cannot
/// be advertised without an answer or answered without being advertised.
const TABLE: [(&str, Evaluate); 4] = [
    ("so3_exp", so3_exp),
    ("so3_log", so3_log),
    ("sen3_exp_n1", sen3_exp_n1),
    ("sen3_log_n1", sen3_log_n1),
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

/// A helicoid quaternion `[w, x, y, z]`, as `tf_tree_math`'s `Quat`: the same order.
pub(crate) fn helicoid_to_tf_tree_quat(q: &[f64]) -> Result<Quat, String> {
    match *q {
        [w, x, y, z] => Ok(Quat::new(w, x, y, z)),
        _ => Err(format!("a quaternion has 4 values, not {}", q.len())),
    }
}

/// `tf_tree_math`'s `Quat` as a helicoid quaternion `[w, x, y, z]`.
pub(crate) fn tf_tree_to_helicoid_quat(q: Quat) -> Vec<f64> {
    vec![q.w, q.x, q.y, q.z]
}

/// A helicoid 3-vector `[x, y, z]`, as `tf_tree_math`'s `Vec3`: the same order.
pub(crate) fn helicoid_to_tf_tree_vec3(v: &[f64]) -> Result<Vec3, String> {
    match *v {
        [x, y, z] => Ok(Vec3::new(x, y, z)),
        _ => Err(format!("a vector has 3 values, not {}", v.len())),
    }
}

/// `tf_tree_math`'s `Vec3` as a helicoid 3-vector `[x, y, z]`.
pub(crate) fn tf_tree_to_helicoid_vec3(v: Vec3) -> Vec<f64> {
    vec![v.x, v.y, v.z]
}

/// A helicoid SE(3) tangent `[φ; ρ]` as `tf_tree_math`'s `[ω, v]` array: the same order.
pub(crate) fn helicoid_to_tf_tree_tangent(tau: &[f64]) -> Result<[f64; 6], String> {
    <[f64; 6]>::try_from(tau).map_err(|_| format!("a tangent has 6 values, not {}", tau.len()))
}

/// `tf_tree_math`'s `[ω, v]` array as a helicoid SE(3) tangent `[φ; ρ]`.
pub(crate) fn tf_tree_to_helicoid_tangent(xi: [f64; 6]) -> Vec<f64> {
    xi.to_vec()
}

/// A helicoid `(q, x)` as `tf_tree_math`'s `Iso3 { q, t }`.
pub(crate) fn helicoid_to_tf_tree_iso3(q: &[f64], x: &[f64]) -> Result<Iso3, String> {
    Ok(Iso3::new(
        helicoid_to_tf_tree_quat(q)?,
        helicoid_to_tf_tree_vec3(x)?,
    ))
}

/// `tf_tree_math`'s `Iso3` as a helicoid `(q, x)`.
pub(crate) fn tf_tree_to_helicoid_iso3(t: Iso3) -> (Vec<f64>, Vec<f64>) {
    (tf_tree_to_helicoid_quat(t.q), tf_tree_to_helicoid_vec3(t.t))
}

fn fields<const N: usize>(entries: [(&str, Vec<f64>); N]) -> Fields {
    entries
        .into_iter()
        .map(|(k, v)| (k.to_string(), v))
        .collect()
}

/// `fn_id` evaluated at `input` by `tf_tree_math`; `None` for an id this runner does not answer.
pub(crate) fn answer(fn_id: &str, input: &Fields) -> Option<Result<Fields, String>> {
    let (_, evaluate) = TABLE.iter().find(|(id, _)| *id == fn_id)?;
    Some(evaluate(input))
}

fn so3_exp(input: &Fields) -> Result<Fields, String> {
    let phi = helicoid_to_tf_tree_vec3(take(input, "phi", 3)?)?;
    Ok(fields([("q", tf_tree_to_helicoid_quat(exp_so3(phi)))]))
}

fn so3_log(input: &Fields) -> Result<Fields, String> {
    let q = helicoid_to_tf_tree_quat(take(input, "q", 4)?)?;
    Ok(fields([("phi", tf_tree_to_helicoid_vec3(log_so3(q)))]))
}

fn sen3_exp_n1(input: &Fields) -> Result<Fields, String> {
    let xi = helicoid_to_tf_tree_tangent(take(input, "tau", 6)?)?;
    let (q, x) = tf_tree_to_helicoid_iso3(exp_se3(xi));
    Ok(fields([("q", q), ("x", x)]))
}

fn sen3_log_n1(input: &Fields) -> Result<Fields, String> {
    let t = helicoid_to_tf_tree_iso3(take(input, "q", 4)?, take(input, "x", 3)?)?;
    Ok(fields([("tau", tf_tree_to_helicoid_tangent(log_se3(t)))]))
}

#[cfg(test)]
mod tests {
    use core::f64::consts::{FRAC_1_SQRT_2, FRAC_2_PI, FRAC_PI_2};

    use tf_tree_math::twist::Twist;

    use super::*;

    /// Far above the rounding of the inputs and of `tf_tree_math`, far below any convention slip.
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

    /// `x` of `Exp([skew_phi(); (1, 2, 3)])`, from `mpmath.expm` of the 4×4 matrix at 60 digits.
    const SKEW_X: [f64; 3] = [0.7494107237383275, 1.85168153974187, 3.157688988882956];

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

    #[test]
    fn a_quaternion_is_w_first_on_both_sides() -> Result<(), String> {
        let q = helicoid_to_tf_tree_quat(&[1.0, 2.0, 3.0, 4.0])?;
        assert_eq!(
            (q.w.to_bits(), q.x.to_bits(), q.y.to_bits(), q.z.to_bits()),
            (
                1.0f64.to_bits(),
                2.0f64.to_bits(),
                3.0f64.to_bits(),
                4.0f64.to_bits()
            )
        );
        assert_eq!(tf_tree_to_helicoid_quat(q), [1.0, 2.0, 3.0, 4.0]);
        assert!(helicoid_to_tf_tree_quat(&[1.0; 3]).is_err());
        Ok(())
    }

    #[test]
    fn a_vector_is_x_y_z_on_both_sides() -> Result<(), String> {
        let v = helicoid_to_tf_tree_vec3(&[1.0, 2.0, 3.0])?;
        assert_eq!(v, Vec3::new(1.0, 2.0, 3.0));
        assert_eq!(
            tf_tree_to_helicoid_vec3(Vec3::new(4.0, 5.0, 6.0)),
            [4.0, 5.0, 6.0]
        );
        assert!(helicoid_to_tf_tree_vec3(&[1.0; 2]).is_err());
        Ok(())
    }

    #[test]
    fn an_isometry_is_q_then_t_on_both_sides() -> Result<(), String> {
        let t = helicoid_to_tf_tree_iso3(&[1.0, 2.0, 3.0, 4.0], &[5.0, 6.0, 7.0])?;
        assert_eq!(t.q, Quat::new(1.0, 2.0, 3.0, 4.0));
        assert_eq!(t.t, Vec3::new(5.0, 6.0, 7.0));
        let (q, x) = tf_tree_to_helicoid_iso3(t);
        assert_eq!((q, x), (vec![1.0, 2.0, 3.0, 4.0], vec![5.0, 6.0, 7.0]));
        assert!(helicoid_to_tf_tree_iso3(&[1.0; 4], &[1.0; 4]).is_err());
        assert!(helicoid_to_tf_tree_iso3(&[1.0; 3], &[1.0; 3]).is_err());
        Ok(())
    }

    #[test]
    fn a_tangent_is_rotation_first_on_both_sides() -> Result<(), String> {
        let xi = helicoid_to_tf_tree_tangent(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0])?;
        assert_eq!(Twist::from_se3(xi).omega, Vec3::new(1.0, 2.0, 3.0));
        assert_eq!(Twist::from_se3(xi).v, Vec3::new(4.0, 5.0, 6.0));
        assert_eq!(
            tf_tree_to_helicoid_tangent(xi),
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        );
        assert!(helicoid_to_tf_tree_tangent(&[0.0; 7]).is_err());
        Ok(())
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
        let out = run("so3_exp", &[("phi", &skew_phi())])?;
        near(&out["q"], &skew_q());
        Ok(())
    }

    /// `Log` reads `q` and `-q` as one rotation with `θ ∈ [0, π]`: the flip on `w < 0`.
    #[test]
    fn so3_log_about_a_skew_axis_and_of_its_negative() -> Result<(), String> {
        let q = skew_q();
        let negative = q.map(|c| -c);
        for q in [q, negative] {
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

    /// `τ = [φ; ρ] = [(π/2)(2, 3, 6)/7; (1, 2, 3)]`: the rotation is `skew_q()` and `x = V ρ`
    /// is `SKEW_X`. Read translation-first, the same six numbers would rotate by `(1, 2, 3)`.
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

    /// `q:nonunit` reaches `tf_tree_math` as the corpus has it: the answer is the library's on
    /// the raw input, not on its normalization. `q` is off the unit sphere by far more than the
    /// corpus's `2⁻⁴⁵`, so that the two differ in the last bits (asserted).
    #[test]
    fn a_non_unit_quaternion_is_not_normalized() -> Result<(), String> {
        let q = [0.50048828125, 0.3, 0.40048828125, 0.69951171875];
        let raw = Quat::new(q[0], q[1], q[2], q[3]);
        assert_ne!(
            bits(&tf_tree_to_helicoid_vec3(log_so3(raw))),
            bits(&tf_tree_to_helicoid_vec3(log_so3(raw.normalize())))
        );
        let out = run("so3_log", &[("q", &q)])?;
        assert_eq!(
            bits(&out["phi"]),
            bits(&tf_tree_to_helicoid_vec3(log_so3(raw)))
        );
        let x = [1.0, 2.0, 3.0];
        let out = run("sen3_log_n1", &[("q", &q), ("x", &x)])?;
        let t = Iso3::new(raw, Vec3::new(1.0, 2.0, 3.0));
        assert_eq!(bits(&out["tau"]), bits(&log_se3(t)));
        assert_ne!(bits(&log_se3(t)), bits(&log_se3(t.normalized())));
        Ok(())
    }

    #[test]
    fn a_missing_or_short_input_is_an_error_and_an_unknown_id_is_none() {
        assert!(run("so3_exp", &[]).is_err());
        assert!(run("so3_exp", &[("phi", &[0.0; 4])]).is_err());
        assert!(answer("so3_jr", &Fields::new()).is_none());
        assert!(answer("sen3_exp_n2", &Fields::new()).is_none());
    }

    /// The harness holds the same list (`RUNNERS` in `xtask`) and fails on a listed id that got no
    /// answer file, so dropping an id here is an error there, not a run with fewer rows.
    #[test]
    fn the_answered_ids_are_the_four_of_the_status_table() {
        let ids: Vec<&str> = supported().collect();
        assert_eq!(ids, ["so3_exp", "so3_log", "sen3_exp_n1", "sen3_log_n1"]);
    }
}
