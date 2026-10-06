//! The exact arms, the series arms and the call-site groups (`super`).
//!
//! An exact arm is the "Exact arm computes" column of `NUMERICS.md` §4 at `z > 0`, in the operand
//! order `docs/maths/coefficients.md` CO.6 measured: `θ²` inside a cancelling sum is `fl(θ̂·θ̂)`,
//! not `z`, which raises the constants. A group evaluates it at the safe argument (`0003` item 3).
//! `z = θ² >= 0` is the domain of every function here, `debug_assert!`ed; `θ²` overflowing is
//! CO.15(c), open. A lane evaluates the series arm it does not select at the raw `z`, so that arm is
//! not finite once its Horner form overflows, far beyond any rotation.

use helicoid_linalg::{Mask, Precision, Real};

use super::generated::{
    A_F32, A_F64, B_F32, B_F64, COS_HALF_F32, COS_HALF_F64, C_F32, C_F64, D_F32, D_F64, E_F32,
    E_F64, K_F32, K_F64, R_F32, R_F64,
};
#[cfg(any(test, feature = "__sweep"))]
use super::generated::{
    SWEPT_A_F32, SWEPT_A_F64, SWEPT_B_F32, SWEPT_B_F64, SWEPT_COS_HALF_F32, SWEPT_COS_HALF_F64,
    SWEPT_C_F32, SWEPT_C_F64, SWEPT_D_F32, SWEPT_D_F64, SWEPT_E_F32, SWEPT_E_F64, SWEPT_K_F32,
    SWEPT_K_F64, SWEPT_R_F32, SWEPT_R_F64,
};
use super::Switch;

/// A series at the width of its table; both widths are exact in the binary64 `Real::lit` takes.
#[derive(Clone, Copy)]
enum Terms<'a> {
    F64(&'a [f64]),
    F32(&'a [f32]),
}

/// One coefficient's `Switch` at one precision, its length erased.
#[derive(Clone, Copy)]
struct Arm<'a> {
    below: f64,
    short_below: f64,
    short_terms: usize,
    terms: Terms<'a>,
}

impl<const M: usize> Switch<f64, M> {
    fn arm(&self) -> Arm<'_> {
        Arm {
            below: self.below,
            short_below: self.short_below,
            short_terms: self.short_terms,
            terms: Terms::F64(&self.series),
        }
    }
}

impl<const M: usize> Switch<f32, M> {
    fn arm(&self) -> Arm<'_> {
        Arm {
            below: f64::from(self.below),
            short_below: f64::from(self.short_below),
            short_terms: self.short_terms,
            terms: Terms::F32(&self.series),
        }
    }
}

impl Arm<'_> {
    fn below<S: Real>(&self) -> S {
        S::lit(self.below)
    }

    fn short_below<S: Real>(&self) -> S {
        S::lit(self.short_below)
    }

    fn horner<S: Real>(&self, z: S) -> S {
        match self.terms {
            Terms::F64(t) => horner(t, z),
            Terms::F32(t) => horner(t, z),
        }
    }

    /// The first `short_terms` terms alone, the second arm of `0039` items 5, 6 and 10. Equal to
    /// [`Self::horner`] at every grid point and every corpus record below `short_below`, by the
    /// sweep's own choice of the prefix, so taking it moves no bit of any measured row (`0047`).
    fn short<S: Real>(&self, z: S) -> S {
        let head = self.short_terms;
        match self.terms {
            Terms::F64(t) => horner(t.get(..head).unwrap_or(t), z),
            Terms::F32(t) => horner(t.get(..head).unwrap_or(t), z),
        }
    }
}

/// `Σ terms[j] zʲ` as `p_j = t_j + z p_{j+1}` from the last term (CO.9); no `mul_add`. Under `Dual`
/// it is exactly the derivative of the polynomial.
fn horner<S: Real, T: Copy + Into<f64>>(terms: &[T], z: S) -> S {
    terms
        .iter()
        .rev()
        .fold(S::zero(), |p, &t| S::lit(t.into()) + z * p)
}

/// The table of `S`'s precision, chosen at monomorphization.
fn table<S: Real, T>(wide: T, narrow: T) -> T {
    match S::PRECISION {
        Precision::F64 => wide,
        Precision::F32 => narrow,
    }
}

fn nonnegative<S: Real>(z: S) {
    // Written as "no lane is negative", not "every lane is `>= 0`": the two differ only on NaN,
    // which this assert is not here to reject. `θ²` and `n²` are dot products, so a *negative*
    // branch variable is a sign error upstream and worth a panic; a NaN one is a NaN input, and
    // every arm below returns NaN for it, which is the answer a value function owes its caller.
    debug_assert!(
        !z.lt(S::zero()).any(),
        "coeffs: the branch variable is negative"
    );
}

/// One group: every member on the short prefix below the group's second switch, else on its series
/// arm below its own first switch, else on its exact arm.
///
/// The outer branch is the **group's** second switch, its members' smallest, so everything under it
/// is the one-arm kernel unchanged (`0047` item 7). The inner branch is on "every member is on its
/// series arm", where no exact arm runs. Off it `exact` runs once for the group, so what its
/// members share (`θ`, a `sin_cos`) is formed once, at `select(all, 1, z)`: `z` is then at or above
/// the group's smallest first switch, which is positive, so a member on its series arm is finite
/// there too. Each member then selects by its own mask; a scalar mask runs the Horner of the
/// members that are small and no other.
fn grouped<S: Real, const G: usize>(
    arms: [Arm<'_>; G],
    exact: impl FnOnce(S) -> [S; G],
    z: S,
) -> [S; G] {
    nonnegative(z);
    // A mask per member instead of this one measured **1.15x to 1.17x slower** than the one-arm
    // kernel at near-identity `θ` and 1.20x to 1.28x slower above the second switch, where no
    // member takes a short arm at all: one comparison and one branch each, on every call
    // (`0047` item 7). It costs no accuracy — between the group's switch and a member's own the
    // member takes the *whole* arm, which is the arm the sweep admitted the prefix to agree with
    // there — so neither reading moves a grid point or a corpus record. A member with no second
    // arm has `0`, which no `z` is below, and takes the group's switch down with it.
    let short_below = arms.iter().fold(f64::INFINITY, |m, a| m.min(a.short_below));
    S::branch(
        z.lt(S::lit(short_below)),
        || arms.map(|a| a.short(z)),
        || {
            let small = arms.map(|a| z.lt(a.below()));
            let all = small.iter().fold(S::zero().le(S::zero()), |m, &s| m.and(s));
            S::branch(
                all,
                || arms.map(|a| a.horner(z)),
                || {
                    let x = exact(S::select(all, S::one(), z));
                    core::array::from_fn(|i| S::branch(small[i], || arms[i].horner(z), || x[i]))
                },
            )
        },
    )
}

/// `(k, cos(θ/2))` at `θ² = z`: `Exp`'s quaternion is `(cos(θ/2), k φ)` (`NUMERICS.md` §3.1).
pub(crate) fn exp_coeffs<S: Real>(z: S) -> (S, S) {
    let arms = table::<S, _>(
        [K_F64.arm(), COS_HALF_F64.arm()],
        [K_F32.arm(), COS_HALF_F32.arm()],
    );
    let [k, cos_half] = grouped(arms, exact_k_cos_half, z);
    (k, cos_half)
}

/// `(a, b)` at `θ² = z`: `J = I ∓ aW + bW²` (`NUMERICS.md` §3.5).
pub(crate) fn jr_coeffs<S: Real>(z: S) -> (S, S) {
    let arms = table::<S, _>([A_F64.arm(), B_F64.arm()], [A_F32.arm(), B_F32.arm()]);
    let [a, b] = grouped(arms, exact_a_b, z);
    (a, b)
}

/// `(a, b)` sharing one `θ = sqrt z`, which is what `grouped`'s contract asks of an exact arm and
/// what `exact_b_d_e` and `exact_k_cos_half` already do. `exact_a(z)` and `exact_b(z)` each take
/// their own root, so the pair cost two; sharing it is a common subexpression on the same `z` and
/// changes no bit. `libm::sqrt` is one instruction only where the target maps it to hardware, so
/// `just no-std` (`thumbv7em-none-eabihf`) and `just wasm` were paying two software roots per
/// `SO3::jr` call.
///
/// The two `sin_cos` are **not** collapsed: `sin θ = 2 sin(θ/2) cos(θ/2)` is a different rounding
/// and so owes a measurement (`0006`), where sharing the root owes none.
fn exact_a_b<S: Real>(z: S) -> [S; 2] {
    let th = z.sqrt();
    let k = half_angle(th).0 / th;
    [S::lit(2.0) * k * k, b_from(th, th.sin_cos().0)]
}

/// `c` at `θ² = z`: `J⁻¹ = I ± W/2 + cW²` (`NUMERICS.md` §3.5).
pub(crate) fn jr_inv_coeff<S: Real>(z: S) -> S {
    nonnegative(z);
    let c = table::<S, _>(C_F64.arm(), C_F32.arm());
    S::branch(
        z.lt(c.short_below()),
        || c.short(z),
        || {
            let small = z.lt(c.below());
            S::branch(
                small,
                || c.horner(z),
                || exact_c(S::select(small, S::one(), z)),
            )
        },
    )
}

/// `(b, d, e)` at `θ² = z`: Barfoot's `Q` block (`NUMERICS.md` §5.3).
// SE_N(3)'s alone: `SEn3::jr` and `jr_inv` are its only consumers (`PHASE3.md` §5).
pub(crate) fn q_coeffs<S: Real>(z: S) -> (S, S, S) {
    let arms = table::<S, _>(
        [B_F64.arm(), D_F64.arm(), E_F64.arm()],
        [B_F32.arm(), D_F32.arm(), E_F32.arm()],
    );
    let [b, d, e] = grouped(arms, exact_b_d_e, z);
    (b, d, e)
}

/// Whether [`log_ratio`] would take its **short** arm at `(n2, w)`: `w > 0` and
/// `s = n²/w² < short_below`, `0047`'s second, shortest series arm.
///
/// Exported because `SO3::geodesic` dispatches its two arms on it (`0051`). It is a predicate
/// [`log_ratio`] already computes and not a second number: `0004` forbids typing a switch point,
/// and a caller that wanted one would have to, so it reads this instead. One definition, so the
/// two cannot drift apart — which is the whole reason this is a function and not a copy of the
/// four lines.
///
/// **`below` is the wrong one of `r`'s two switches for that caller**, which is measured and not
/// argued: dispatching on it leaves 50 of `so3_geodesic/geo:generic`'s 60 records on the provided
/// body, where it reads `2.5019 u` against oracle #1's `1.834` and loses the stratum
/// (`measure_geodesic::scan`'s table, `0051`). `short_below` is four decades of `s` lower and
/// dominates every stratum, with every smaller threshold reading identically — so the choice is
/// not delicate, but it is a choice, and
/// `measure_geodesic::tests::the_shipped_dispatch_threshold_still_dominates` fails if a later sweep
/// moves this number into the region that loses.
#[inline]
pub(crate) fn log_ratio_takes_short_arm<S: Real>(n2: S, w: S) -> S::Mask {
    nonnegative(n2);
    let r = table::<S, _>(R_F64.arm(), R_F32.arm());
    let positive = S::zero().lt(w);
    // `n² < T w²`, not `n²/w² < T`: the same test for `w > 0` and a multiply instead of a divide,
    // which is 1.1 ns [`log_ratio`] would then spend again on its own copy of the quotient. The two
    // can part by an ulp exactly at the boundary, and nothing reads them as one number: a caller
    // dispatching on this picks an arm, both arms are correct, and `0051` measured five decades of
    // slack either side. `w * w` cannot overflow for a quaternion within a factor of two of unit,
    // and where it does the comparison is false, which is the arm that assumes nothing.
    positive.and(n2.lt(r.short_below::<S>() * w * w))
}

/// `r = 2 atan2(n, w)/n` at `n² = n2` (`NUMERICS.md` §3.2). The series arm is taken iff `w > 0` and
/// `s = n²/w² < switch`: the reading in `super`. Each arm is at its safe argument, the series arm's
/// `w` (`CO.16(d)`) too, and the mask forms `s` from `w² = 1` where `w` is not positive, so no
/// operation is non-finite there. Where `w > 0` and `s` overflows, `super` says what a lane sees.
pub(crate) fn log_ratio<S: Real>(n2: S, w: S) -> S {
    nonnegative(n2);
    let r = table::<S, _>(R_F64.arm(), R_F32.arm());
    let positive = S::zero().lt(w);
    let s = n2 / S::select(positive, w * w, S::one());
    // The second switch is the **outermost** branch, so everything under it is what it was
    // (`0047` item 7): measured the other way round — the prefix selected inside the series arm —
    // this function was **1.16x slower** at `θ = 0.5`, where `s = 0.065` takes the whole arm, even
    // though that spelling also saved a division. One comparison on the way in buys the prefix;
    // a comparison *under* the series arm buys it and pays for two Horner bodies in one block.
    let tiny = positive.and(s.lt(r.short_below()));
    S::branch(
        tiny,
        || {
            // `tiny` implies `w > 0`, so the selects are the safe arguments a lane that is not
            // selected needs (`0003` item 3): `2 / 1` and a Horner at `0`, both finite. `s` is
            // the Horner's argument here because it is `n2 / (w * w)` to the bit on a selected
            // lane — the mask's own `select` took `w * w` — and a division is 1.11 ns.
            S::lit(2.0) / S::select(tiny, w, S::one()) * r.short(S::select(tiny, s, S::zero()))
        },
        || {
            let small = positive.and(s.lt(r.below()));
            S::branch(
                small,
                || {
                    let w = S::select(small, w, S::one());
                    S::lit(2.0) / w * r.horner(n2 / (w * w))
                },
                || exact_r(S::select(small, S::one(), n2), w),
            )
        },
    )
}

/// `sin(θ/2)` and `cos(θ/2)` from one `sin_cos`.
fn half_angle<S: Real>(th: S) -> (S, S) {
    (S::lit(0.5) * th).sin_cos()
}

// The per-coefficient arms the sweep and the tests scan one at a time; the shipped groups take
// `exact_a_b`, `exact_k_cos_half` and `exact_b_d_e`, which share what their members share.
#[cfg(any(test, feature = "__sweep"))]
pub(super) fn exact_k<S: Real>(z: S) -> S {
    let th = z.sqrt();
    half_angle(th).0 / th
}

/// `(k, cos(θ/2))` from the one `sin_cos` of `θ/2` that `NUMERICS.md` §3.1 states.
fn exact_k_cos_half<S: Real>(z: S) -> [S; 2] {
    let th = z.sqrt();
    let (s, c) = half_angle(th);
    [s / th, c]
}

/// `2k²`, exact where `(1 - cos θ)/θ²` is not (`NUMERICS.md` §4).
#[cfg(any(test, feature = "__sweep"))]
pub(super) fn exact_a<S: Real>(z: S) -> S {
    let k = exact_k(z);
    S::lit(2.0) * k * k
}

#[cfg(any(test, feature = "__sweep"))]
pub(super) fn exact_b<S: Real>(z: S) -> S {
    let th = z.sqrt();
    b_from(th, th.sin_cos().0)
}

/// `(θ - sin θ)/θ³` from `θ` and `s = sin θ`.
fn b_from<S: Real>(th: S, s: S) -> S {
    let t2 = th * th;
    (th - s) / (t2 * th)
}

/// `1/θ² - cot(θ/2)/(2θ)`, regular at `π`.
pub(super) fn exact_c<S: Real>(z: S) -> S {
    let th = z.sqrt();
    let t2 = th * th;
    let (s, co) = half_angle(th);
    S::one() / t2 - co / (S::lit(2.0) * th * s)
}

/// `(θ² - 4 sin²(θ/2))/(2θ⁴)` from `θ`.
fn d_from<S: Real>(th: S) -> S {
    let t2 = th * th;
    let s = half_angle(th).0;
    (t2 - S::lit(4.0) * s * s) / (S::lit(2.0) * t2 * t2)
}

/// `(2θ - 3 sin θ + θ cos θ)/(2θ⁵)` from `θ`, `s = sin θ` and `co = cos θ`.
fn e_from<S: Real>(th: S, s: S, co: S) -> S {
    let t2 = th * th;
    ((S::lit(2.0) * th - S::lit(3.0) * s) + th * co) / (S::lit(2.0) * th * (t2 * t2))
}

/// `(b, d, e)`, `sin θ` and `cos θ` formed once for `b` and `e`.
fn exact_b_d_e<S: Real>(z: S) -> [S; 3] {
    let th = z.sqrt();
    let (s, co) = th.sin_cos();
    [b_from(th, s), d_from(th), e_from(th, s, co)]
}

#[cfg(any(test, feature = "__sweep"))]
pub(super) fn exact_d<S: Real>(z: S) -> S {
    d_from(z.sqrt())
}

#[cfg(any(test, feature = "__sweep"))]
pub(super) fn exact_e<S: Real>(z: S) -> S {
    let th = z.sqrt();
    let (s, co) = th.sin_cos();
    e_from(th, s, co)
}

/// `cos(θ/2)`, the quaternion's scalar part.
#[cfg(any(test, feature = "__sweep"))]
pub(super) fn exact_cos_half<S: Real>(z: S) -> S {
    half_angle(z.sqrt()).1
}

/// `2 atan2(n, w)/n` at `n² = n2 > 0`.
pub(super) fn exact_r<S: Real>(n2: S, w: S) -> S {
    let n = n2.sqrt();
    S::lit(2.0) * n.atan2(w) / n
}

/// The first `terms` of a swept series: how many the sweep asks for is its own, at most all.
#[cfg(any(test, feature = "__sweep"))]
fn swept<S: Real>(wide: &[f64], narrow: &[f32], z: S, terms: usize) -> S {
    debug_assert!((1..=wide.len()).contains(&terms), "coeffs: 1..=8 terms");
    match S::PRECISION {
        Precision::F64 => horner(wide.get(..terms).unwrap_or(wide), z),
        Precision::F32 => horner(narrow.get(..terms).unwrap_or(narrow), z),
    }
}

#[cfg(any(test, feature = "__sweep"))]
pub(super) fn series_k<S: Real>(z: S, terms: usize) -> S {
    swept(&SWEPT_K_F64, &SWEPT_K_F32, z, terms)
}

#[cfg(any(test, feature = "__sweep"))]
pub(super) fn series_a<S: Real>(z: S, terms: usize) -> S {
    swept(&SWEPT_A_F64, &SWEPT_A_F32, z, terms)
}

#[cfg(any(test, feature = "__sweep"))]
pub(super) fn series_b<S: Real>(z: S, terms: usize) -> S {
    swept(&SWEPT_B_F64, &SWEPT_B_F32, z, terms)
}

#[cfg(any(test, feature = "__sweep"))]
pub(super) fn series_c<S: Real>(z: S, terms: usize) -> S {
    swept(&SWEPT_C_F64, &SWEPT_C_F32, z, terms)
}

#[cfg(any(test, feature = "__sweep"))]
pub(super) fn series_d<S: Real>(z: S, terms: usize) -> S {
    swept(&SWEPT_D_F64, &SWEPT_D_F32, z, terms)
}

#[cfg(any(test, feature = "__sweep"))]
pub(super) fn series_e<S: Real>(z: S, terms: usize) -> S {
    swept(&SWEPT_E_F64, &SWEPT_E_F32, z, terms)
}

#[cfg(any(test, feature = "__sweep"))]
pub(super) fn series_cos_half<S: Real>(z: S, terms: usize) -> S {
    swept(&SWEPT_COS_HALF_F64, &SWEPT_COS_HALF_F32, z, terms)
}

/// `2/w` times the series in `s = n²/w²`, at `w > 0`.
#[cfg(any(test, feature = "__sweep"))]
pub(super) fn series_r<S: Real>(n2: S, w: S, terms: usize) -> S {
    S::lit(2.0) / w * swept(&SWEPT_R_F64, &SWEPT_R_F32, n2 / (w * w), terms)
}
