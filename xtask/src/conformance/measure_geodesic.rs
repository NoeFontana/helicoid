//! A measurement, not a routine ([`0049`] *Further work* 1): which spelling of `SO3::geodesic` is
//! the most accurate, per stratum, over the committed `so3_geodesic` corpus.
//!
//! `helicoid` won `geo:consecutive` against `tf_tree_math::slerp` and lost `geo:generic` and
//! `geo:near-pi`, and `0006` reads the maximum per stratum, so the question was whether the losses
//! were the *route* or the arithmetic. Seven spellings of the one function GE.14 proves them all
//! equal to, scored by the same exact-rational metric the conformance runner uses
//! ([`metric::rule`]), so the figures are comparable with the committed rows without a conversion:
//!
//! - [`route_s`] is what **ships**, and it is the **control**: its figures have to reproduce the
//!   committed `helicoid` rows, or the harness is wrong and nothing else here means anything.
//! - [`route_a`] is the **provided body**, `q0 Exp(t Log(q0* q1))`, which shipped alone until
//!   `0050` and was `0051`'s arm below the switch until `0059` wrote it out as [`route_f`].
//! - [`route_b`] is GE.14's grouped middle expression, `q0 (cos t·alpha, varpi_t v)`. This is the
//!   rotation part of GE.12, i.e. of `PHASE4.md` §1.2's screw twin, so what it measures is whether
//!   that twin can fix these two strata.
//! - [`route_c`] is GE.14's right-hand expression, the blend `[sin((1-t)a) q0 + sin(ta) q1]/sin a`,
//!   with `alpha` from the quaternion product so no `sqrt(1 - d^2)` appears. This is the spelling
//!   `tf_tree_math::slerp` computes, which is the oracle that wins these two strata.
//! - [`route_f`] is the provided body as `0059` writes it out, the arm below the switch that
//!   ships: the same function, with `r` and `Exp`'s coefficients no longer waiting on each other's
//!   products.
//! - [`route_d`] is the same blend with `alpha` from the 4-dot alone,
//!   `atan2(sqrt(1 - d^2), d)`, which is the cheap form: one dot instead of a Hamilton product.
//!   It is here to price the cancellation in `1 - d^2`, which is total where the two rotations are
//!   consecutive.
//!
//! No `acos` anywhere, in any route: `Real` has no `acos` (`0022`) and D5 bans it on this path, so
//! a candidate spelled with one could not ship and measuring it would answer nothing.
//!
//! [`scan`] then answers the question `0051` needed: at which threshold a *two-arm* routine takes
//! the provided body, measured per stratum rather than assumed.
//!
//! Scored at **binary64 only**, which is where the three committed domination failures are, and
//! that is this module's standing gap: `0051` ships a switch whose `f32` boundary is
//! `R_F32.short_below`, 29x the `f64` one, and **no scan covers it**, because `so3_geodesic` has no
//! `f32` stratum to score. [`THRESHOLDS`] carries the `f32` value so the gap is at least bracketed
//! at `f64`; closing it needs an `@f32` geodesic stratum (`0016`), which `0051` *Further work* 2
//! owes. `0006` reads the bars per precision and this module reads one.
//!
//! [`se3`] is the same question for `SE3::geodesic` over `se3_geodesic` (`0050` *Further work* 2),
//! asked of `PHASE4.md` §1.2's screw twin (`0054`).
//!
//! [`0049`]: ../../../docs/decisions/0049-the-boundary-is-what-removes-a-way-to-be-wrong.md

use std::collections::BTreeMap;

use helicoid::{LieGroup, Quat, SEn3, SE3, SO3};
use helicoid_linalg::Precision;

use super::corpus::{parse_line, Record};
use super::metric::{self, Score};
use super::subject::Output;

/// `q` as the corpus holds it, `w` first, in the shape [`metric`] scores.
fn out_of(q: [f64; 4]) -> Output {
    Output::from([("q".to_string(), q.to_vec())])
}

/// Input `key` as a quaternion, `w` first; nothing when the record does not hold four values.
fn quat_of(rec: &Record, key: &str) -> Option<[f64; 4]> {
    rec.input(key).and_then(|v| v.first_chunk::<4>()).copied()
}

/// The Hamilton product, `w` first, as `Quat`'s `Mul` writes it.
fn mul(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    let ([aw, ax, ay, az], [bw, bx, by, bz]) = (a, b);
    [
        aw * bw - ax * bx - ay * by - az * bz,
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
    ]
}

/// `q0* q1` with `Log`'s sign rule applied (`w >= 0`, GE.14's hypothesis), and the flipped `q1`.
///
/// The flip is on the *product*, not on `q0 . q1`: they are the same number — `w(q0* q1)` **is**
/// the 4-dot — but taking it from the product keeps one spelling for every route below, so a sign
/// difference between two of them cannot be the flip.
fn rel(q0: [f64; 4], q1: [f64; 4]) -> ([f64; 4], [f64; 4]) {
    let conj = [q0[0], -q0[1], -q0[2], -q0[3]];
    let d = mul(conj, q1);
    // `copysign`, not `d[0] < 0.0`: the shipped `geodesic` and `SO3::log` both flip on
    // `1.0.copysign(w)`, which is `-1` at `w = -0.0` where `<` is false. A `-0.0` product is
    // reachable, and the two spellings would then blend `q1` against `-q1` -- a different
    // quaternion, and a bit-equality test comparing two different functions.
    let flip = 1.0_f64.copysign(d[0]);
    (d.map(|x| flip * x), q1.map(|x| flip * x))
}

/// Lift a stored quaternion into `SO3` as the subject does.
///
/// `from_quat_unchecked` for `shipped.rs`'s stated reason: the reference is the geodesic between the
/// two rotations the records *denote*, and the half `u` of `|‖q‖² − 1|` a program carrying the
/// stored quaternion reads for it is the shipped behaviour the row measures.
fn lift(q: [f64; 4]) -> SO3<f64> {
    let [w, x, y, z] = q;
    SO3::from_quat_unchecked(Quat { w, x, y, z })
}

/// `w` first, as the corpus holds it.
fn wxyz(q: &Quat<f64>) -> [f64; 4] {
    [q.w, q.x, q.y, q.z]
}

/// The **provided body**, `q0 Exp(t Log(q0* q1))`, which is what `SO3::geodesic` was until `0050`.
///
/// Taken through `reference::geodesic`, which *is* that body — `x0.rplus(&x1.rminus(x0).scale(t))`
/// — so this column keeps measuring the route the override replaced, and keeps doing so after the
/// override. Reading it from `SO3::geodesic` instead is what [`route_s`] is for, and conflating the
/// two is what made both of this module's tests fire the moment the override landed.
fn route_a(q0: [f64; 4], q1: [f64; 4], t: f64) -> [f64; 4] {
    wxyz(&helicoid::reference::geodesic(&lift(q0), &lift(q1), t).quat())
}

/// What **ships**: `SO3::geodesic` itself, which `0050` made [`route_e`].
///
/// Its column is the control, and it being a separate column from [`route_a`] is the point: the one
/// asserts the harness against the committed rows, the other keeps the superseded route measurable.
fn route_s(q0: [f64; 4], q1: [f64; 4], t: f64) -> [f64; 4] {
    wxyz(&SO3::geodesic(&lift(q0), &lift(q1), t).quat())
}

/// GE.14's grouped form: `q0 (cos t·alpha, varpi_t v)` with `varpi_t = sin(t·alpha)/sin(alpha)`.
///
/// `sin(alpha)` is `‖v‖` exactly for a unit quaternion, so the division is by the norm that was
/// measured and not by a `sin` recomputed from `alpha` — one rounding fewer and no second
/// transcendental. This is the rotation part of the screw form (GE.12), which is the whole point of
/// measuring it: it keeps both Hamilton products that [`route_a`] has and replaces only the
/// coefficient evaluation.
fn route_b(q0: [f64; 4], q1: [f64; 4], t: f64) -> [f64; 4] {
    let (d, _) = rel(q0, q1);
    let v = [d[1], d[2], d[3]];
    let nv = libm::sqrt(v.iter().map(|x| x * x).sum::<f64>());
    let alpha = libm::atan2(nv, d[0]);
    let (s, c) = libm::sincos(t * alpha);
    // `0/0` at two bitwise equal rotations: the arc is a point and the answer is `q0` itself.
    let w = if nv == 0.0 { 0.0 } else { s / nv };
    mul(q0, [c, w * v[0], w * v[1], w * v[2]])
}

/// GE.14's right-hand form: `[sin((1-t)alpha) q0 + sin(t alpha) q1] / sin(alpha)`.
///
/// `alpha` comes from the quaternion product, as [`route_b`]'s does, so the chord `1 - d^2` never
/// appears. One Hamilton product instead of [`route_a`]'s two, and a linear blend of two unit
/// quaternions instead of the second one: every coefficient is `O(1)` and both weights are
/// non-negative on `[0, 1]`, so the blend cannot cancel there.
fn route_c(q0: [f64; 4], q1: [f64; 4], t: f64) -> [f64; 4] {
    let (d, q1) = rel(q0, q1);
    let nv = libm::sqrt(d[1..].iter().map(|x| x * x).sum::<f64>());
    let alpha = libm::atan2(nv, d[0]);
    blend(q0, q1, t, alpha, nv)
}

/// [`route_c`]'s blend with `alpha` from the 4-dot alone: `atan2(sqrt(1 - d^2), d)`.
///
/// The cheap form — one dot, no Hamilton product — and the one that pays for it: `1 - d^2` loses
/// every digit of `‖v‖` as `d` approaches 1, which is the whole of `geo:consecutive`.
fn route_d(q0: [f64; 4], q1: [f64; 4], t: f64) -> [f64; 4] {
    let (d, q1) = rel(q0, q1);
    let dot = d[0];
    let nv = libm::sqrt((1.0 - dot * dot).max(0.0));
    let alpha = libm::atan2(nv, dot);
    blend(q0, q1, t, alpha, nv)
}

/// [`route_c`]'s blend dividing by `sin(alpha)` **recomputed**, not by `‖v‖`.
///
/// One change, and it is the one the oracle makes: every weight then comes from the same `alpha`
/// through the same `sin`, so `sin(t a)/sin(a)` is a number divided by **itself** at `t = 1` and the
/// numerator is an exact zero at `t = 0`. Both endpoints are reproduced bit for bit. Dividing by
/// `‖v‖` instead mixes a `sqrt`-derived denominator with `sin`-derived numerators, so the ratio at
/// `t = 1` is `1 + eps` and the endpoint is off by `eps` -- which is why [`route_c`]'s and
/// [`route_d`]'s worst `geo:generic` records are endpoints.
fn route_e(q0: [f64; 4], q1: [f64; 4], t: f64) -> [f64; 4] {
    let (d, q1) = rel(q0, q1);
    let nv = libm::sqrt(d[1..].iter().map(|x| x * x).sum::<f64>());
    let alpha = libm::atan2(nv, d[0]);
    blend(q0, q1, t, alpha, libm::sincos(alpha).0)
}

/// [`route_e`] with `alpha` from the **chord**, so no Hamilton product and no `1 - d^2`.
///
/// `h = ½‖q0 - q1‖²` is `1 - |q0 . q1|` computed from component differences, which is the trick
/// `tf_tree_math::slerp` opens with and the one this measurement was missing: each difference is
/// exact by Sterbenz while the two quaternions are within a factor of two of each other, so `h`
/// carries no cancellation however consecutive they are -- where `1 - d^2` has lost every digit
/// ([`route_d`] reads `5.1e7 u` at `geo:consecutive` for exactly that). Then `d = 1 - h` and
/// `‖v‖ = sin(alpha) = sqrt(h(2 - h))`.
///
/// Shippable as written: `atan2`, `sqrt`, `sin`, no `acos` (`0022`, D5).
fn route_g(q0: [f64; 4], q1: [f64; 4], t: f64) -> [f64; 4] {
    let dot = q0.iter().zip(q1).map(|(a, b)| a * b).sum::<f64>();
    let q1 = if dot < 0.0 { q1.map(|x| -x) } else { q1 };
    let h = 0.5
        * q0.iter()
            .zip(q1)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>();
    let nv = libm::sqrt(h * (2.0 - h));
    let alpha = libm::atan2(nv, 1.0 - h);
    blend(q0, q1, t, alpha, libm::sincos(alpha).0)
}

/// The provided body as `0059` writes it out, which is `SO3::geodesic` below `r`'s second switch:
/// `θ² = r² (t² n²)` and `q0 (c, k t r v) = c q0 + k t r (q0 (0, v))`, with `r` and `(k, c)` the
/// kernel's own (`__sweep`). A transcription, so [`tests::the_shipped_route_is_the_arm_the_rule_picks`]
/// ties the shipped arm to a stated expression rather than to itself.
fn route_f(q0: [f64; 4], q1: [f64; 4], t: f64) -> [f64; 4] {
    let (d, _) = rel(q0, q1);
    let n2 = (d[1] * d[1] + d[2] * d[2]) + d[3] * d[3];
    let r = helicoid::__sweep::log_ratio(n2, d[0]);
    let p = mul(q0, [0.0, d[1], d[2], d[3]]);
    let (k, c) = helicoid::__sweep::exp_coeffs(r * r * (t * t * n2));
    let tr = t * r;
    [0, 1, 2, 3].map(|i| c * q0[i] + k * (tr * p[i]))
}

/// The shared tail of the four blends: the two sines over `den`, which is `sin(alpha)` by one
/// spelling or another and is what the routes above differ in.
///
/// At `den == 0` the two rotations are equal and the blend is `0/0`; the arc is a point, so the
/// answer is `q0`, which is what both endpoints of a degenerate geodesic are.
fn blend(q0: [f64; 4], q1: [f64; 4], t: f64, alpha: f64, den: f64) -> [f64; 4] {
    if den == 0.0 {
        return q0;
    }
    // `sincos(..).0` and not `sin(..)`: the shipped routine goes through `Real::sin_cos`, and a
    // bit-equality test between the two would be asserting that `libm`'s `sin` and `sincos` agree
    // -- true today, promised by nothing. Same call, no assumption.
    let (s0, s1) = (
        libm::sincos((1.0 - t) * alpha).0 / den,
        libm::sincos(t * alpha).0 / den,
    );
    [0, 1, 2, 3].map(|i| s0 * q0[i] + s1 * q1[i])
}

/// The corpus this measurement reads, relative to `xtask`'s manifest.
const CORPUS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../conformance/corpus/so3_geodesic.jsonl"
);

/// The committed rows this measurement is read against (`PHASE4.md` §0.0), in `u`, so the control
/// is checked by the test and not by a reader.
///
/// `shipped` is [`route_s`]'s row. It read 1.572 / 2.721 / 2.429 while the provided body shipped
/// alone — [`route_a`]'s column — and 1.644 / 1.738 / 1.642 under `0050`'s single blend; `0051`'s
/// two arms took the better cell of each, and `0059`'s written-out arm reads 0.993 at
/// `geo:consecutive`, where it runs.
/// The oracle column is at **full precision**, not the three decimals the prose quotes: the shipped
/// route reads `1.6417` at `geo:near-pi` against an oracle of `1.6417180974225531`, so a bound of
/// `1.642` would leave less margin than its own rounding and could not tell dominating from losing.
const COMMITTED: [(&str, f64, f64); 3] = [
    // stratum, shipped (= `0051`'s two arms), `tf_tree_math::slerp`
    ("geo:consecutive", 0.993, 2.1873313039389336),
    ("geo:generic", 1.738, 1.8341289492373978),
    ("geo:near-pi", 1.642, 1.6417180974225531),
];

/// One spelling: `(q0, q1, t)` to the quaternion at parameter `t`.
type Route = fn([f64; 4], [f64; 4], f64) -> [f64; 4];

/// The spellings scored, in the order every table below prints them. [`route_s`] is last so the
/// six studied spellings keep their columns and the shipped one is read beside them.
const ROUTES: [Route; 8] = [
    route_a, route_b, route_c, route_d, route_e, route_f, route_g, route_s,
];

/// [`route_s`]'s index in [`ROUTES`]: the control's column, which is the last by that array's own
/// stated invariant — a literal `6` would silently repoint both controls at `route_g` the moment a
/// route is appended, which is the failure the controls exist to prevent.
const SHIPPED: usize = ROUTES.len() - 1;

/// One route's worst reading on one stratum: the maximum in `u`, and the record that attained it.
type Worst = (f64, u64);

/// Every route's [`Worst`] per stratum, in [`ROUTES`]' order.
type Table = BTreeMap<String, Vec<Worst>>;
/// `R_F64`'s **second** switch, which `SO3::geodesic` dispatches on (`0051`).
///
/// A literal here and generated there, which is the only way a test outside the crate can state the
/// rule: `coeffs` is `pub(crate)`, so this is a transcription, and
/// [`tests::the_shipped_route_is_the_arm_the_rule_picks`] fails the moment the two disagree — which
/// is exactly what a sweep moving the number should do.
const R_SHORT_F64: f64 = 5.048065716667471e-5;

/// Candidate dispatch thresholds on `s = tan²α = n²/w²`, for [`scan`].
///
/// The two `R_F64` arms are the only *generated* numbers a two-arm `geodesic` could borrow without
/// a sweep of its own (`0051`); the decades are there to show where the accuracy crossover is, so
/// the two can be compared rather than one assumed.
const THRESHOLDS: [(&str, f64); 9] = [
    ("R_F64.below (θ≈5.8e-1)", 8.976871324473142e-2),
    // `R_F32.short_below`, the boundary the `f32` build actually ships, 29x the `f64` one. Scored
    // here at `f64` because no `f32` stratum exists to score it at: it brackets the gap rather than
    // closing it, and it is the only entry between the two above and below it.
    ("R_F32.short (θ≈7.7e-2)", 1.485508e-3),
    ("R_F64.short (θ≈1.4e-2)", R_SHORT_F64),
    ("s<1e-6  (θ≈2.0e-3)", 1e-6),
    ("s<1e-7  (θ≈6.3e-4)", 1e-7),
    ("s<1e-8  (θ≈2.0e-4)", 1e-8),
    ("s<1e-10 (θ≈2.0e-5)", 1e-10),
    ("s<1e-12 (θ≈2.0e-6)", 1e-12),
    ("s<0     (blend only)", 0.0),
];

/// Every route's maximum `u` per stratum over the committed corpus, plus the worst record's id.
fn measure() -> Result<Table, String> {
    let text = std::fs::read_to_string(CORPUS).map_err(|e| format!("{CORPUS}: {e}"))?;
    let rule = metric::rule("so3_geodesic").ok_or("no metric rule for `so3_geodesic`")?;
    let mut out: BTreeMap<String, Vec<(f64, u64)>> = BTreeMap::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let rec = parse_line(line)?;
        let (Some(q0), Some(q1)) = (quat_of(&rec, "q0"), quat_of(&rec, "q1")) else {
            return Err(format!("record {} holds no `q0`/`q1`", rec.id));
        };
        let &[t] = rec
            .input("t")
            .and_then(|v| v.first_chunk::<1>())
            .ok_or_else(|| format!("record {} holds no `t`", rec.id))?;
        let row = out
            .entry(rec.stratum.clone())
            .or_insert_with(|| vec![(0.0, 0); ROUTES.len()]);
        for (slot, route) in row.iter_mut().zip(ROUTES) {
            let got = out_of(route(q0, q1, t));
            match rule.score(&rec, &got, Precision::F64)? {
                Score::Finite(u) if u > slot.0 => *slot = (u, rec.id),
                Score::Finite(_) => {}
                other => return Err(format!("record {} scored {other:?}", rec.id)),
            }
        }
    }
    Ok(out)
}

/// For each [`THRESHOLDS`] entry, the per-stratum max of a two-arm routine that takes [`route_f`]
/// where `s < T` and [`route_e`] elsewhere: [`route_a`] until `0059`, the same function.
///
/// This is the evidence a dispatch threshold needs and the thing `0050` *Further work* 1 asserted
/// without: that the provided body is both the faster *and* the more accurate arm below some
/// switch. It is faster up to `R_F64.below`; whether it is more accurate that far is this scan's
/// question, and the answer decides whether a two-arm needs a sweep of its own.
fn scan() -> Result<BTreeMap<String, Vec<f64>>, String> {
    let text = std::fs::read_to_string(CORPUS).map_err(|e| format!("{CORPUS}: {e}"))?;
    let rule = metric::rule("so3_geodesic").ok_or("no metric rule for `so3_geodesic`")?;
    let mut out: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let rec = parse_line(line)?;
        let (Some(q0), Some(q1)) = (quat_of(&rec, "q0"), quat_of(&rec, "q1")) else {
            return Err(format!("record {} holds no `q0`/`q1`", rec.id));
        };
        let &[t] = rec
            .input("t")
            .and_then(|v| v.first_chunk::<1>())
            .ok_or_else(|| format!("record {} holds no `t`", rec.id))?;
        // `s` from the same product the shipped predicate reads. Not the same *expression*:
        // `coeffs::log_ratio_takes_short_arm` spells the test `n2 < T w²` to avoid forming the
        // quotient at all, and the two can part by an ulp exactly at the boundary. That is the
        // point of the scan rather than a flaw in it -- it asks which decade the threshold belongs
        // in, and the answer has five decades of slack.
        let (d, _) = rel(q0, q1);
        let n2 = d[1..].iter().map(|x| x * x).sum::<f64>();
        let s = n2 / (d[0] * d[0]);
        let row = out
            .entry(rec.stratum.clone())
            .or_insert_with(|| vec![0.0; THRESHOLDS.len()]);
        for (slot, (_, limit)) in row.iter_mut().zip(THRESHOLDS) {
            // `t >= 1` goes to the blend on every threshold, as the shipped dispatch does: only
            // the blend is exact there, and that is not a property to trade for a threshold.
            let q = if s < limit && t < 1.0 {
                route_f(q0, q1, t)
            } else {
                route_e(q0, q1, t)
            };
            match rule.score(&rec, &out_of(q), Precision::F64)? {
                Score::Finite(u) if u > *slot => *slot = u,
                Score::Finite(_) => {}
                other => return Err(format!("record {} scored {other:?}", rec.id)),
            }
        }
    }
    Ok(out)
}

/// Which spelling of `SE3::geodesic` is the most accurate, per stratum, over `se3_geodesic`
/// (`0054`), by the metric the conformance runner uses.
///
/// - [`se3::route_a`], the **provided body**, `X₀ Exp(t Log(X₀⁻¹X₁))`, which shipped until `0054`.
/// - [`se3::route_b`], GE.12 **verbatim**: `tf_tree_math::dualquat::screw_pow`'s exact arm, one
///   `atan2` and one `sin_cos`, `ϖ_t = sin(tα)/‖v‖`, the rotation `q₀ q_rᵗ`.
/// - [`se3::route_c`], `route_b`'s translation under the shipped rotation, `SO3::geodesic`'s:
///   between `b` and `c` only the rotation moves, so the pair prices it.
/// - [`se3::route_s`], what **ships**: `route_c`'s rotation, the definition's translation below
///   `r`'s second switch and GE.12's above it, both in the world frame. The control, which has
///   to reproduce the subject's rows.
pub(super) mod se3 {
    use super::*;

    /// The corpus this half reads, relative to `xtask`'s manifest.
    pub(super) const CORPUS: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../conformance/corpus/se3_geodesic.jsonl"
    );

    /// One spelling: `(X₀, X₁, t)` to `(q, x)`, `w` first.
    pub(super) type Route = fn(&SE3<f64>, &SE3<f64>, f64) -> ([f64; 4], [f64; 3]);

    /// The spellings, in the order the table prints them; [`route_s`] last, as above.
    pub(super) const ROUTES: [Route; 4] = [route_a, route_b, route_c, route_s];

    /// `X` as the subject reads it, `from_quat_unchecked` for `lift`'s reason.
    pub(super) fn pose(rec: &Record, qk: &str, xk: &str) -> Option<SE3<f64>> {
        let q = quat_of(rec, qk)?;
        let x = *rec.input(xk)?.first_chunk::<3>()?;
        Some(SEn3::from_parts(lift(q), [helicoid_linalg::Vector(x)]))
    }

    fn split(x: &SE3<f64>) -> ([f64; 4], [f64; 3]) {
        (wxyz(&x.rotation().quat()), x.translation().0)
    }

    /// The provided body, through `reference::geodesic`, for [`route_a`]'s reason in the SO(3) half.
    pub(super) fn route_a(x0: &SE3<f64>, x1: &SE3<f64>, t: f64) -> ([f64; 4], [f64; 3]) {
        split(&helicoid::reference::geodesic(x0, x1, t))
    }

    /// What ships.
    pub(super) fn route_s(x0: &SE3<f64>, x1: &SE3<f64>, t: f64) -> ([f64; 4], [f64; 3]) {
        split(&SE3::geodesic(x0, x1, t))
    }

    /// GE.12 as `screw_pow`'s exact arm writes it, with the rotation `q₀ q_rᵗ`.
    pub(super) fn route_b(x0: &SE3<f64>, x1: &SE3<f64>, t: f64) -> ([f64; 4], [f64; 3]) {
        let (qrt, xt) = power(x0, x1, t);
        let q0 = wxyz(&x0.rotation().quat());
        (mul(q0, qrt), xt)
    }

    /// [`route_b`]'s translation, the shipped rotation.
    pub(super) fn route_c(x0: &SE3<f64>, x1: &SE3<f64>, t: f64) -> ([f64; 4], [f64; 3]) {
        let (_, xt) = power(x0, x1, t);
        let q = SO3::geodesic(&x0.rotation(), &x1.rotation(), t).quat();
        (wxyz(&q), xt)
    }

    /// GE.12 verbatim: `(q_rᵗ, x₀ + R₀ x_Δᵗ)`, `Δ` through `inv_mul` as every route forms it.
    ///
    /// The `0/0` at `‖v‖ = 0` is `screw_pow`'s degenerate arm, the provided body; no corpus
    /// record reaches it, so it is there for the arithmetic's sake and not measured.
    fn power(x0: &SE3<f64>, x1: &SE3<f64>, t: f64) -> ([f64; 4], [f64; 3]) {
        let delta = x1.inv_mul(x0);
        let dq = wxyz(&delta.rotation().quat());
        let flip = 1.0_f64.copysign(dq[0]);
        let [w, x, y, z] = dq.map(|c| flip * c);
        let n2 = (x * x + y * y) + z * z;
        if n2 == 0.0 {
            let (q, xt) = route_a(x0, x1, t);
            let q0 = wxyz(&x0.rotation().quat());
            return (mul([q0[0], -q0[1], -q0[2], -q0[3]], q), xt);
        }
        let sh = libm::sqrt(n2);
        let alpha = libm::atan2(sh, w);
        let (s, c) = libm::sincos(t * alpha);
        let varpi = s / sh;
        let [px, py, pz] = delta.translation().0;
        // `q_d = ½ (0, x_Δ) ⊗ (w, v)`.
        let qd = mul([0.0, px, py, pz], [w, x, y, z]).map(|e| 0.5 * e);
        let kappa = qd[0] / n2;
        let m = [
            qd[1] + x * (kappa * w),
            qd[2] + y * (kappa * w),
            qd[3] + z * (kappa * w),
        ];
        let qrt = [c, varpi * x, varpi * y, varpi * z];
        let k2 = (t * kappa) * c;
        let qdt = [
            (t * varpi) * qd[0],
            m[0] * varpi - x * k2,
            m[1] * varpi - y * k2,
            m[2] * varpi - z * k2,
        ];
        let tv = mul(qdt, [qrt[0], -qrt[1], -qrt[2], -qrt[3]]);
        let xt = helicoid_linalg::Vector([2.0 * tv[1], 2.0 * tv[2], 2.0 * tv[3]]);
        let out = x0.translation() + x0.rotation().act(xt);
        (qrt, out.0)
    }

    /// Every route's [`Worst`] per stratum, in [`ROUTES`]' order.
    pub(super) fn measure() -> Result<Table, String> {
        let text = std::fs::read_to_string(CORPUS).map_err(|e| format!("{CORPUS}: {e}"))?;
        let rule = metric::rule("se3_geodesic").ok_or("no metric rule for `se3_geodesic`")?;
        let mut out = Table::new();
        for line in text.lines().filter(|l| !l.trim().is_empty()) {
            let rec = parse_line(line)?;
            let (Some(x0), Some(x1)) = (pose(&rec, "q0", "x0"), pose(&rec, "q1", "x1")) else {
                return Err(format!("record {} holds no `X0`/`X1`", rec.id));
            };
            let &[t] = rec
                .input("t")
                .and_then(|v| v.first_chunk::<1>())
                .ok_or_else(|| format!("record {} holds no `t`", rec.id))?;
            let row = out
                .entry(rec.stratum.clone())
                .or_insert_with(|| vec![(0.0, 0); ROUTES.len()]);
            for (slot, route) in row.iter_mut().zip(ROUTES) {
                let (q, x) = route(&x0, &x1, t);
                let got =
                    Output::from([("q".to_string(), q.to_vec()), ("x".to_string(), x.to_vec())]);
                match rule.score(&rec, &got, Precision::F64)? {
                    Score::Finite(u) if u > slot.0 => *slot = (u, rec.id),
                    Score::Finite(_) => {}
                    other => return Err(format!("record {} scored {other:?}", rec.id)),
                }
            }
        }
        Ok(out)
    }

    /// The subject's rows (`0054`, `geo:consecutive` `0059`'s) and `tf_tree_math`'s `ScLerp`, at
    /// full precision.
    pub(super) const COMMITTED: [(&str, f64, f64); 3] = [
        ("geo:consecutive", 0.993, 2.3363543564610385),
        ("geo:generic", 2.057, 2.501902372122162),
        ("geo:near-pi", 2.721, 3.2527931951103226),
    ];
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Bit equality up to the sign of a zero, which is the exception the shipped `geodesic`'s
    /// rustdoc already names: a `-0.0` component of an endpoint comes back `+0.0` because a weight
    /// of exactly zero times it is added to it.
    #[allow(clippy::float_cmp)]
    fn same(a: &[f64; 4], b: &[f64; 4]) -> bool {
        a.iter()
            .zip(b)
            .all(|(&x, &y)| x.to_bits() == y.to_bits() || (x == 0.0 && y == 0.0))
    }

    /// `(q0, q1, t)` of one corpus line.
    fn case(line: &str) -> Result<([f64; 4], [f64; 4], f64, u64), String> {
        let rec = parse_line(line)?;
        let (Some(q0), Some(q1)) = (quat_of(&rec, "q0"), quat_of(&rec, "q1")) else {
            return Err(format!("record {} holds no `q0`/`q1`", rec.id));
        };
        let &[t] = rec
            .input("t")
            .and_then(|v| v.first_chunk::<1>())
            .ok_or_else(|| format!("record {} holds no `t`", rec.id))?;
        Ok((q0, q1, t, rec.id))
    }

    /// Every non-blank corpus line.
    fn lines() -> Result<Vec<String>, String> {
        let text = std::fs::read_to_string(CORPUS).map_err(|e| format!("{CORPUS}: {e}"))?;
        Ok(text
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(str::to_string)
            .collect())
    }

    /// The control, which runs in `just test`: [`route_s`] -- what ships -- reproduces the
    /// committed `helicoid` row of every stratum to the three decimals `PHASE4.md` §0.0 quotes,
    /// **and** dominates oracle #1 on each.
    ///
    /// Without the first half the figures the measurement prints are unreadable: a harness that
    /// mis-scores would mis-score every route and still look self-consistent. The second half is
    /// the bar `0006` reads, and it is a check rather than a hope because `SO3::geodesic` dispatches
    /// on a number the **sweep** owns: `r`'s second switch is chosen to minimise a corpus-weighted
    /// term count, not to place this boundary, so a later sweep may move it. [`scan`]'s table says
    /// which move would matter -- `R_F64.below`, four decades of `s` higher, reads `2.5019 u` at
    /// `geo:generic` and loses it -- and every threshold at or below the shipped one reads
    /// identically, so there are five decades of slack.
    ///
    /// One `measure()` for both halves: it parses 180 records and scores seven routes through the
    /// exact-rational metric, 1260 scorings, and two tests asking about the same column would run
    /// all of it twice in `just test`.
    #[test]
    fn the_shipped_route_reproduces_the_committed_rows_and_dominates() -> Result<(), String> {
        let got = measure()?;
        for (stratum, shipped, oracle) in COMMITTED {
            let (u, id) = *got
                .get(stratum)
                .and_then(|r| r.get(SHIPPED))
                .ok_or_else(|| format!("no `{stratum}` in the corpus"))?;
            assert!(
                (u - shipped).abs() < 5e-4,
                "{stratum}: the shipped route reads {u} u, `PHASE4.md` §0.0 says {shipped}"
            );
            assert!(
                u <= oracle,
                "{stratum}: the shipped geodesic reads {u} u at record #{id}, over oracle #1's \
                 {oracle}"
            );
        }
        Ok(())
    }

    /// `SO3::geodesic` **is** one of the two routes on every corpus record, bit for bit, and which
    /// one follows `0051`'s rule: [`route_f`] below the switch (`0059`), [`route_e`] above.
    ///
    /// The tables alone cannot say this: columns agreeing to three decimals is not one being the
    /// other. This is what ties the records' measurements to the routine they decided — and it
    /// pins the dispatch rule itself, since it recomputes the predicate from `s` and `t` rather
    /// than trusting whichever arm happens to match.
    #[test]
    fn the_shipped_route_is_the_arm_the_rule_picks() -> Result<(), String> {
        let (mut f_arm, mut e_arm) = (0usize, 0usize);
        for line in lines()? {
            let (q0, q1, t, id) = case(&line)?;
            let (d, _) = rel(q0, q1);
            let n2 = d[1..].iter().map(|x| x * x).sum::<f64>();
            let s = n2 / (d[0] * d[0]);
            // `0051`: the provided body below `r`'s second switch, except at `t >= 1`.
            let provided = d[0] > 0.0 && s < R_SHORT_F64 && t < 1.0;
            let want = if provided {
                f_arm += 1;
                route_f(q0, q1, t)
            } else {
                e_arm += 1;
                route_e(q0, q1, t)
            };
            let got = route_s(q0, q1, t);
            assert!(
                got.iter()
                    .zip(&want)
                    .all(|(a, b)| a.to_bits() == b.to_bits()),
                "record {id} at t = {t}, s = {s:e}: shipped {got:?} is not the \
                 {} arm's {want:?}",
                if provided { "provided" } else { "blend" }
            );
        }
        // Both arms are exercised, or the test is checking one of them against itself.
        assert!(
            f_arm > 0 && e_arm > 0,
            "arms taken: {f_arm} provided, {e_arm} blend"
        );
        Ok(())
    }

    /// Why [`route_e`] wins the two losing strata: it reproduces **both** endpoints bit for bit,
    /// where the **provided** body ([`route_a`]) reproduces only `t = 0`.
    ///
    /// At `t = 0` every route's left weight is `sin(a)/den`; only when `den` is that same `sin(a)`
    /// is the ratio exactly one. At `t = 1` the mirror holds. The corpus pins this rather than the
    /// argument doing: every record at `t = 0` or `t = 1`, all three strata, and the provided body
    /// asserted to **fail** at `t = 1` on at least one of them -- if it ever stops failing, the
    /// asymmetry `PHASE4.md` §0.0 records has gone and this reasoning needs rereading.
    #[test]
    fn the_recomputed_denominator_is_what_makes_both_endpoints_exact() -> Result<(), String> {
        let (mut ends, mut a_fails_at_one) = (0usize, 0usize);
        for line in lines()? {
            let (q0, q1, t, id) = case(&line)?;
            // `rel`'s flip: at `t = 1` the answer is `q1` up to the sign `Log` chose, which is the
            // same rotation and what the oracle documents for its own `s = 1`.
            let want = if t.to_bits() == 0.0f64.to_bits() {
                q0
            } else if t.to_bits() == 1.0f64.to_bits() {
                rel(q0, q1).1
            } else {
                continue;
            };
            ends += 1;
            let got = route_e(q0, q1, t);
            assert!(
                same(&got, &want),
                "route E, record {id} at t = {t}: {got:?} is not {want:?}"
            );
            if t.to_bits() == 1.0f64.to_bits() && !same(&route_a(q0, q1, t), &want) {
                a_fails_at_one += 1;
            }
        }
        assert!(
            ends >= 30,
            "only {ends} endpoint records; the corpus changed"
        );
        assert!(
            a_fails_at_one > 0,
            "the provided body is now exact at t = 1 on every record; \
             `PHASE4.md` §0.0's one-sided claim and this test's reasoning both need rereading"
        );
        Ok(())
    }

    /// `R_SHORT_F64` is still the number the shipped routine switches on, to the bit.
    ///
    /// The transcription is the weak point of every claim in this module: `coeffs` is `pub(crate)`,
    /// so nothing here can read the generated constant, and
    /// [`the_shipped_route_is_the_arm_the_rule_picks`] only notices a drift if some corpus record's
    /// `s` happens to land between the old threshold and the new one. A sweep moving it a few
    /// percent would leave no record in the gap and every test green.
    ///
    /// So this one builds the pair instead of waiting for one: `q0 = 1` and
    /// `q1 = (w, nv, 0, 0)` with `nv = sqrt(s/(1+s))`, `w = sqrt(1/(1+s))` realise any `s` to
    /// rounding. One `s` a hair under the threshold has to take the provided body and one a hair
    /// over has to take the blend; the arms differ by about 7 `u`, so the comparison is decisive,
    /// and the test says so rather than assuming it.
    #[test]
    fn the_shipped_threshold_is_where_the_arms_actually_part() -> Result<(), String> {
        let pair = |s: f64| -> ([f64; 4], [f64; 4]) {
            let (nv, w) = ((s / (1.0 + s)).sqrt(), (1.0 / (1.0 + s)).sqrt());
            ([1.0, 0.0, 0.0, 0.0], [w, nv, 0.0, 0.0])
        };
        // `t` interior, so neither endpoint's exactness can hide which arm ran.
        let t = 1.0 / 3.0;
        for (name, s, want_provided) in [
            ("just under", R_SHORT_F64 * (1.0 - 1e-9), true),
            ("just over", R_SHORT_F64 * (1.0 + 1e-9), false),
        ] {
            let (q0, q1) = pair(s);
            let got = route_s(q0, q1, t);
            let same = |x: &[f64; 4]| got.iter().zip(x).all(|(p, q)| p.to_bits() == q.to_bits());
            let (took_f, took_e) = (same(&route_f(q0, q1, t)), same(&route_e(q0, q1, t)));
            assert!(
                took_f != took_e,
                "{name} the threshold at s = {s:e}: the arms agree to the bit, so this test \
                 cannot see which ran -- choose a `t` or an `s` where they differ"
            );
            assert_eq!(
                took_f,
                want_provided,
                "{name} the threshold at s = {s:e}: the shipped routine took the {} arm, so \
                 `R_SHORT_F64` is no longer `R_F64.short_below` and `scan`'s table is stale",
                if took_f { "provided" } else { "blend" }
            );
        }
        Ok(())
    }

    /// Where a two-arm's threshold may sit: the per-stratum max for each [`THRESHOLDS`] entry.
    ///
    /// `cargo test -p xtask -- --ignored --nocapture scan_the_dispatch_threshold`.
    #[test]
    #[ignore = "a measurement, printed for a record to cite"]
    #[allow(clippy::print_stdout)]
    fn scan_the_dispatch_threshold() -> Result<(), String> {
        let got = scan()?;
        println!(
            "so3_geodesic, binary64: max u per stratum of `route_f where s < T, else route_e`"
        );
        println!(
            "{:<26} {:>16} {:>14} {:>14}   verdict",
            "T", "consecutive", "generic", "near-pi"
        );
        for (i, (name, _)) in THRESHOLDS.iter().enumerate() {
            let cell = |st: &str| got.get(st).and_then(|r| r.get(i)).copied().unwrap_or(0.0);
            let (c, g, n) = (
                cell("geo:consecutive"),
                cell("geo:generic"),
                cell("geo:near-pi"),
            );
            // Domination is per stratum against oracle #1, which is what `0006` reads.
            let beats = c <= 2.187 && g <= 1.834 && n <= 1.642;
            println!(
                "{name:<26} {c:>16.4} {g:>14.4} {n:>14.4}   {}",
                if beats { "dominates" } else { "LOSES" }
            );
        }
        Ok(())
    }

    /// `cargo test -p xtask -- --ignored --nocapture measure_the_geodesic_spellings`.
    #[test]
    #[ignore = "a measurement, printed for a record to cite; the two tests above are the checks"]
    #[allow(clippy::print_stdout)]
    fn measure_the_geodesic_spellings() -> Result<(), String> {
        const NAMES: [&str; ROUTES.len()] = [
            "A provided",
            "B grouped",
            "C blend/|v|",
            "D dot/|v|",
            "E blend/sin",
            "F written out",
            "G chord/sin",
            "S shipped",
        ];
        let got = measure()?;
        println!("so3_geodesic, binary64, max u per stratum (worst record id)");
        print!("{:<16}", "stratum");
        for n in NAMES {
            print!("{n:>21}");
        }
        println!("{:>10}", "oracle");
        for (stratum, _, oracle) in COMMITTED {
            let row = got
                .get(stratum)
                .ok_or_else(|| format!("no `{stratum}` in the corpus"))?;
            print!("{stratum:<16}");
            for &(u, id) in row {
                print!("{u:>13.3} u (#{id:>3})");
            }
            println!("{oracle:>10.3}");
        }
        Ok(())
    }

    /// [`se3::route_s`] reproduces the subject's rows and dominates `tf_tree_math`'s `ScLerp` on
    /// every stratum: the SE(3) half's control, for the SO(3) control's reasons (`0054`).
    #[test]
    fn the_shipped_se3_route_reproduces_its_rows_and_dominates() -> Result<(), String> {
        let got = se3::measure()?;
        let last = se3::ROUTES.len() - 1;
        for (stratum, shipped, oracle) in se3::COMMITTED {
            let (u, id) = *got
                .get(stratum)
                .and_then(|r| r.get(last))
                .ok_or_else(|| format!("no `{stratum}` in the corpus"))?;
            assert!(
                (u - shipped).abs() < 5e-4,
                "{stratum}: the shipped route reads {u} u, `0054` says {shipped}"
            );
            assert!(
                u <= oracle,
                "{stratum}: the shipped geodesic reads {u} u at record #{id}, over `ScLerp`'s {oracle}"
            );
        }
        Ok(())
    }

    /// `cargo test -p xtask -- --ignored --nocapture measure_the_se3_geodesic_spellings`.
    #[test]
    #[ignore = "a measurement, printed for a record to cite"]
    #[allow(clippy::print_stdout)]
    fn measure_the_se3_geodesic_spellings() -> Result<(), String> {
        const NAMES: [&str; se3::ROUTES.len()] =
            ["A provided", "B GE.12", "C GE.12+SO3", "S shipped"];
        let got = se3::measure()?;
        println!("se3_geodesic, binary64, max u per stratum (worst record id)");
        print!("{:<16}", "stratum");
        for n in NAMES {
            print!("{n:>21}");
        }
        println!("{:>10}", "oracle");
        for (stratum, _, oracle) in se3::COMMITTED {
            let row = got
                .get(stratum)
                .ok_or_else(|| format!("no `{stratum}` in the corpus"))?;
            print!("{stratum:<16}");
            for &(u, id) in row {
                print!("{u:>13.3} u (#{id:>3})");
            }
            println!("{oracle:>10.3}");
        }
        Ok(())
    }
}
