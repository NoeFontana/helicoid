//! The Sim(3) coefficient group against `conformance/corpus/sim3_ref.jsonl` and closed forms at
//! large `|σ|` (`0070`): the shipped arms meet the sweep's recorded objective at both precisions,
//! the value-only call is the full one's `v` to the bit, and the scaled arrangement holds to
//! `ln MAX`.

use super::sim3::{sim3_coeffs, sim3_coeffs_scaled, sim3_coeffs_scaled_pair, sim3_v, Sim3Coeffs};
use helicoid_linalg::{Dual, Real};
use std::format;
use std::string::{String, ToString};
use std::vec::Vec;

const REF: &str = include_str!("../../../../conformance/corpus/sim3_ref.jsonl");
const SWEEP: &str = include_str!("../../../../conformance/sweeps/thresholds-sim3.csv");
const NAMES: [&str; 8] = ["v0", "v1", "v2", "w0", "w1", "w2", "dv1", "dv2"];
const WEIGHT: [i32; 8] = [0, 0, 1, 0, 1, 2, 2, 3];

/// `0x1.8p-3`-style hex floats, exactly: the mantissa's digits as an integer, scaled by `2^p`.
fn hex(s: &str) -> Result<f64, String> {
    let (neg, s) = s.strip_prefix('-').map_or((false, s), |r| (true, r));
    let s = s
        .strip_prefix("0x")
        .ok_or_else(|| format!("`{s}` is not hex"))?;
    let (mant, exp) = s
        .split_once('p')
        .ok_or_else(|| format!("`{s}` has no exponent"))?;
    let (int, frac) = mant.split_once('.').unwrap_or((mant, ""));
    let digits = u64::from_str_radix(&[int, frac].concat(), 16).map_err(|e| e.to_string())?;
    let e: i32 = exp.parse().map_err(|e| format!("{e}"))?;
    // At most 14 hex digits, so `digits` is exact in an `f64` and the scaling is a power of two.
    let x = (digits as f64) * 2f64.powi(e - 4 * frac.len() as i32);
    Ok(if neg { -x } else { x })
}

fn field<'a>(line: &'a str, key: &str) -> Result<&'a str, String> {
    let at = line
        .find(&format!("\"{key}\":\""))
        .ok_or_else(|| format!("no {key}"))?;
    let rest = &line[at + key.len() + 4..];
    Ok(&rest[..rest.find('"').ok_or("unterminated")?])
}

struct Point {
    sigma: f64,
    z: f64,
    refs: [f64; 8],
}

fn points(precision: &str) -> Result<Vec<Point>, String> {
    let mut out = Vec::new();
    for line in REF.lines() {
        if field(line, "precision")? != precision {
            continue;
        }
        let mut refs = [0.0; 8];
        for (r, n) in refs.iter_mut().zip(NAMES) {
            *r = field(line, n)?.parse().map_err(|e| format!("{n}: {e}"))?;
        }
        out.push(Point {
            sigma: hex(field(line, "sigma")?)?,
            z: hex(field(line, "z")?)?,
            refs,
        });
    }
    Ok(out)
}

fn recorded(precision: &str) -> Result<f64, String> {
    let mut lines = SWEEP.lines();
    let header: Vec<&str> = lines.next().ok_or("empty")?.split(',').collect();
    let col = header
        .iter()
        .position(|h| *h == "objective")
        .ok_or("no objective")?;
    let row = lines
        .map(|l| l.split(',').collect::<Vec<_>>())
        .find(|c| c[0] == precision)
        .ok_or("no row")?;
    row[col].parse().map_err(|e| format!("{e}"))
}

/// The worst weighted error over the rays, as the sweep scores it but against each reference read
/// once into an `f64`: within the recorded objective plus that rounding, at most `|c| u` and so at
/// most `1` in these units at binary64, and nothing at binary32.
fn worst<S: Real + Into<f64>>(
    precision: &str,
    to: impl Fn(f64) -> S,
    u: f64,
) -> Result<f64, String> {
    let mut worst = 0.0f64;
    for p in points(precision)? {
        let c = sim3_coeffs(to(p.sigma), to(p.z));
        let got = [
            c.v[0], c.v[1], c.v[2], c.w[0], c.w[1], c.w[2], c.dv[0], c.dv[1],
        ]
        .map(Into::into);
        let (theta, scale) = (p.z.sqrt(), p.refs[0].abs().max(1.0));
        for j in 0..8 {
            let e = (got[j] - p.refs[j]).abs() * theta.powi(WEIGHT[j]) / (u * scale);
            worst = if e.is_nan() {
                f64::INFINITY
            } else {
                worst.max(e)
            };
        }
        let (v, chi_delta) = sim3_v(to(p.sigma), to(p.z));
        let (v, chi_delta): ([f64; 3], f64) = (v.map(Into::into), chi_delta.into());
        if v.map(f64::to_bits) != [got[0], got[1], got[2]].map(f64::to_bits)
            || chi_delta.to_bits() != c.chi_delta.into().to_bits()
        {
            return Err(format!(
                "sim3_v is not sim3_coeffs' v and chi_delta at ({}, {})",
                p.sigma, p.z
            ));
        }
    }
    Ok(worst)
}

#[test]
fn the_rays_meet_the_recorded_objective() -> Result<(), String> {
    let w64 = worst::<f64>("f64", |x| x, f64::EPSILON / 2.0)?;
    let w32 = worst::<f32>("f32", |x| x as f32, f64::from(f32::EPSILON) / 2.0)?;
    let (r64, r32) = (recorded("f64")?, recorded("f32")?);
    assert!(
        w64 <= r64 + 1.0 && w32 <= r32 + 1e-6,
        "f64 {w64} vs {r64}, f32 {w32} vs {r32}"
    );
    Ok(())
}

/// At `Δ = 0` the series is the limits, exactly.
#[test]
fn the_origin_is_the_limits() {
    let c = sim3_coeffs(0.0_f64, 0.0);
    let want = [
        1.0,
        0.5,
        1.0 / 6.0,
        0.5,
        1.0 / 6.0,
        1.0 / 24.0,
        -1.0 / 24.0,
        -1.0 / 120.0,
    ];
    let got = [
        c.v[0], c.v[1], c.v[2], c.w[0], c.w[1], c.w[2], c.dv[0], c.dv[1],
    ];
    assert_eq!(got.map(f64::to_bits), want.map(f64::to_bits));
}

/// `(σ, z, [v₀ … ∂_z v₂])` from the closed forms at `400 + 2|σ|` digits (mpmath 1.3.0), past the
/// rays: the scaled arrangement for `σ > 0`, the plain one for `σ < 0` (SM.7(g)).
const LARGE: [(f64, f64, [&str; 8]); 12] = [
    (
        10.0,
        1.0,
        [
            "2.2025465794806716517e+3",
            "1.7172907760680670328e+3",
            "8.4082247596788737301e+2",
            "2.2015465794806716517e+2",
            "1.5655627383334540643e+2",
            "6.8426620213454196658e+1",
            "-2.5364735514494786132e+2",
            "-6.0455224850526559817e+1",
        ],
    ),
    (
        -10.0,
        1.0,
        [
            "9.9995460007023751515e-2",
            "9.8969647829424614184e-3",
            "9.8760945399168987169e-4",
            "9.0000453999297624849e-2",
            "7.9408168476984487542e-3",
            "6.9532073937067588825e-4",
            "-9.7123744918674367764e-5",
            "-9.5354863696835395717e-6",
        ],
    ),
    (
        100.0,
        1.0,
        [
            "2.6881171418161354484e+41",
            "2.2472238973116110188e+41",
            "1.2132490126760882134e+41",
            "2.6881171418161354484e+39",
            "2.2322519908211284337e+39",
            "1.1909264927678769291e+39",
            "-3.9366383733102318452e+40",
            "-1.0079832397438834616e+40",
        ],
    ),
    (
        100.0,
        9.0,
        [
            "2.6881171418161354484e+41",
            "1.5292356436886487626e+40",
            "5.9284109126710524907e+40",
            "2.6881171418161354484e+39",
            "1.7923677806448495773e+38",
            "5.9104872348646039949e+38",
            "-1.5411461279228092155e+40",
            "-5.7475054774777248414e+39",
        ],
    ),
    (
        -100.0,
        1.0,
        [
            "1.0e-2",
            "9.9990000999900009999e-5",
            "9.9990000999900009999e-7",
            "9.9e-3",
            "9.7990400939908008999e-5",
            "9.6990500929909008899e-7",
            "-9.9980002999600049994e-9",
            "-9.9980002999600049994e-11",
        ],
    ),
    (
        300.0,
        1.0,
        [
            "6.4747546508041864553e+127",
            "5.4365966828793638486e+127",
            "2.9583077941012915242e+127",
            "2.1582515502680621518e+125",
            "1.8082716392119046007e+125",
            "9.8007502590305749274e+124",
            "-9.6596579571005948205e+126",
            "-2.4905081085766912295e+126",
        ],
    ),
    (
        -300.0,
        9.0,
        [
            "3.3333333333333333333e-3",
            "1.1110000111100001111e-5",
            "3.703333370366667037e-8",
            "3.3222222222222222222e-3",
            "1.1035940849618815889e-5",
            "3.6663025053050497395e-8",
            "-1.2343210246864203703e-10",
            "-4.1144034156214012343e-13",
        ],
    ),
    (
        700.0,
        1e-8,
        [
            "1.4489029353357207278e+301",
            "1.4468330715950083962e+301",
            "7.2238456267615793857e+300",
            "2.069861336193886754e+298",
            "2.0639474432415782913e+298",
            "1.0290294503327376575e+298",
            "-2.4045184436864694158e+300",
            "-6.0027452963928707136e+299",
        ],
    ),
    (
        700.0,
        1.0,
        [
            "1.4489029353357207278e+301",
            "1.21808894313942572e+301",
            "6.6431721133308693388e+300",
            "2.069861336193886754e+298",
            "1.7385223182679124536e+298",
            "9.4654098430688431632e+297",
            "-2.1731326974433946424e+300",
            "-5.6142000922508030121e+299",
        ],
    ),
    (
        700.0,
        9.0,
        [
            "1.4489029353357207278e+301",
            "7.0204255696761754202e+299",
            "3.2026703816397466928e+300",
            "2.069861336193886754e+298",
            "1.0321540909630900501e+297",
            "4.5737688964982622896e+297",
            "-8.3425451842564361144e+299",
            "-3.1690724224460216298e+299",
        ],
    ),
    (
        -700.0,
        1.0,
        [
            "1.4285714285714285714e-3",
            "2.0408121616078334534e-6",
            "2.9154459451540477906e-9",
            "1.4265306122448979592e-3",
            "2.0349812816172804406e-6",
            "2.9029511938173234183e-9",
            "-4.1649142789664377285e-12",
            "-5.9498775413806253264e-15",
        ],
    ),
    (
        -700.0,
        9.0,
        [
            "1.4285714285714285714e-3",
            "2.0407788428375805342e-6",
            "2.9153983469108293346e-9",
            "1.4265306122448979592e-3",
            "2.0349481532380576423e-6",
            "2.9029039355587811614e-9",
            "-4.1647782853734942302e-12",
            "-5.9496832648192774717e-15",
        ],
    ),
];

/// Each scalar's relative error in `u` at the `LARGE` points, finite where the unscaled forms
/// overflow from `σ ≈ 703` (SM.7(g)).
fn large_errors() -> Result<Vec<(f64, f64, [f64; 8])>, String> {
    let mut out = Vec::new();
    for (sigma, z, refs) in LARGE {
        let c = sim3_coeffs(sigma, z);
        let got = [
            c.v[0], c.v[1], c.v[2], c.w[0], c.w[1], c.w[2], c.dv[0], c.dv[1],
        ];
        let mut e = [0.0f64; 8];
        for ((x, g), r) in e.iter_mut().zip(got).zip(refs) {
            let r: f64 = r.parse().map_err(|e| format!("{e}"))?;
            let err = (g - r).abs() / (r.abs() * f64::EPSILON / 2.0);
            *x = if err.is_finite() { err } else { f64::INFINITY };
        }
        out.push((sigma, z, e));
    }
    Ok(out)
}

/// Six scalars within `5u` at every `|σ|`; `w₁ = v₂ - σw₂` and `∂_z v₁ = e^σ a' - σ ∂_z v₂` subtract
/// `σ` times a term of their own size and grow with `σ > 0`, measured at most `2.32σ u` (`w₁` at
/// `σ = 700`, `θ = 3`), against which the `Q` block's own assembly loses `0.6σ² u` (SM.7(g)(iii)):
/// bounded at `3σ u` (`0070`).
#[test]
fn the_scaled_arrangement_holds_to_ln_max() -> Result<(), String> {
    for (sigma, z, e) in large_errors()? {
        for (j, &x) in e.iter().enumerate() {
            let bound = if j == 4 || j == 6 {
                (3.0 * sigma).max(5.0)
            } else {
                5.0
            };
            assert!(x <= bound, "σ = {sigma}, z = {z}, scalar {j}: {x} u");
        }
    }
    Ok(())
}

#[test]
#[ignore = "measurement: prints the figures the large-σ bounds are recorded from"]
#[allow(clippy::print_stdout)]
fn measure_sim3_large() -> Result<(), String> {
    for (sigma, z, e) in large_errors()? {
        std::println!("σ = {sigma}, z = {z}: {e:.1?} u");
    }
    let w = worst::<f64>("f64", |x| x, f64::EPSILON / 2.0)?;
    std::println!("rays f64 {w:.3} u");
    Ok(())
}

/// Under `Dual` the value is the plain one's, bit for bit, on both arms.
#[test]
fn dual_value_is_plain_value() -> Result<(), String> {
    for p in points("f64")?.iter().step_by(37) {
        let plain = sim3_coeffs(p.sigma, p.z);
        let d = sim3_coeffs(Dual::<f64, 2>::variable(p.sigma, 0), Dual::variable(p.z, 1));
        let a = [plain.v, plain.w].concat();
        let b: Vec<f64> = [d.v, d.w].concat().iter().map(|x| x.v).collect();
        if a.iter()
            .map(|x| x.to_bits())
            .ne(b.iter().map(|x| x.to_bits()))
        {
            return Err(String::from("dual value differs"));
        }
    }
    Ok(())
}

fn bits<S: Real + Into<f64>>(c: &Sim3Coeffs<S>) -> [u64; 8] {
    let x = [
        c.v[0], c.v[1], c.v[2], c.w[0], c.w[1], c.w[2], c.dv[0], c.dv[1],
    ];
    x.map(|x| Into::<f64>::into(x).to_bits())
}

/// The pair against two calls at `(σ, z)`; `Err` names the first difference.
fn pair_is_two_calls<S: Real + Into<f64>>(sigma: S, z: S) -> Result<(), String> {
    let pair = sim3_coeffs_scaled_pair(sigma, z);
    for ((c, scaled), s) in pair.iter().zip([sigma, -sigma]) {
        let (want, want_scaled) = sim3_coeffs_scaled(s, z);
        let flag = |m| Into::<f64>::into(S::select(m, S::one(), S::zero())).to_bits();
        if bits(c) != bits(&want) || flag(*scaled) != flag(want_scaled) {
            let (s, z): (f64, f64) = (s.into(), z.into());
            return Err(format!("σ = {s:e}, z = {z:e}"));
        }
    }
    Ok(())
}

/// `sim3_coeffs_scaled_pair` is two `sim3_coeffs_scaled` calls to the bit (`0073`), on both arms,
/// at `σ = ±0`, across `form_b`, `σ_w` and the scaled arrangement, at both precisions and on
/// `Dual`'s value and derivatives: the parity argument holds for any `Real` whose negation is
/// exact.
#[test]
fn the_pair_is_two_calls_to_the_bit() -> Result<(), String> {
    let sigmas = [
        0.0, -0.0, 1e-300, 1e-9, 1e-4, 0.01, 0.1, 0.3, 0.5, 1.0, 2.0, 3.0, 10.0, 80.0, 700.0,
    ];
    let zs = [0.0, 1e-12, 1e-6, 0.01, 0.1, 0.3, 1.0, 2.5, 4.0, 9.8];
    for &s in &sigmas {
        for &z in &zs {
            for s in [s, -s] {
                pair_is_two_calls(s, z)?;
                pair_is_two_calls(s as f32, z as f32)?;
                let (ds, dz) = (Dual::<f64, 2>::variable(s, 0), Dual::variable(z, 1));
                let pair = sim3_coeffs_scaled_pair(ds, dz);
                for ((c, _), s) in pair.iter().zip([ds, -ds]) {
                    let (want, _) = sim3_coeffs_scaled(s, dz);
                    let flat = |c: &Sim3Coeffs<Dual<f64, 2>>| {
                        [c.v, c.w]
                            .concat()
                            .into_iter()
                            .chain(c.dv)
                            .flat_map(|x| [x.v, x.d[0], x.d[1]].map(f64::to_bits))
                            .collect::<Vec<_>>()
                    };
                    if flat(c) != flat(&want) {
                        return Err(format!("Dual at σ = {:e}, z = {z:e}", s.v));
                    }
                }
            }
        }
    }
    for precision in ["f64", "f32"] {
        for p in points(precision)? {
            match precision {
                "f64" => pair_is_two_calls(p.sigma, p.z)?,
                _ => pair_is_two_calls(p.sigma as f32, p.z as f32)?,
            }
        }
    }
    Ok(())
}

/// `sinc_a_b_d`'s members are their own groups' to the bit, on both arms and across every switch:
/// `alpha_value`'s `sin θ/θ`, `jr_coeffs`' `a, b`, `gamma2_coeffs`' `d` (grouping is a cost, never a
/// value).
#[test]
fn sinc_a_b_d_is_its_members_own_groups() {
    fn case<S: Real>(z: S) {
        let (s, a, b, d) = super::kernel::sinc_a_b_d(z);
        let (ja, jb) = super::jr_coeffs(z);
        let (_, gd) = super::gamma2_coeffs(z);
        let want = [super::kernel::alpha_value(z), ja, jb, gd];
        let got = [s, a, b, d];
        assert!(
            got.iter()
                .zip(&want)
                .all(|(g, w)| g.value_f64().to_bits() == w.value_f64().to_bits()),
            "z = {}",
            z.value_f64()
        );
    }
    let mut z = 1e-16_f64;
    while z < 9.8 {
        case(z);
        case(z as f32);
        z *= 1.0371;
    }
}

/// `χΔ` (`0096`, draft): at `σ = 0` near `π`, where `v₀ − zv₂` cancels to `sin θ/θ`, `χ` is
/// `sin θ/θ` to a few `u` of itself; away from `π` it is `v₀ − zv₂` to a few `u` of `v₀`. The
/// subtraction reads about `10⁶ u` of `χ` at `π − 10⁻⁶`.
#[test]
fn chi_delta_is_v0_minus_z_v2_without_the_cancellation() {
    let u = f64::EPSILON / 2.0;
    for theta in [
        2.9,
        3.0,
        3.1,
        core::f64::consts::PI - 1e-6,
        core::f64::consts::PI - 1e-9,
    ] {
        let z = theta * theta;
        let (_, chi_delta) = sim3_v(0.0, z);
        let (chi, sinc) = (chi_delta / z, theta.sin() / theta);
        assert!(
            (chi - sinc).abs() <= 4.0 * u * sinc,
            "theta {theta}: {chi} against {sinc}"
        );
    }
    for (sigma, theta) in [
        (0.5, 1.0),
        (-0.5, 1.0),
        (3.0, 1.5),
        (-3.0, 0.7),
        (1e-3, 2.0),
    ] {
        let z: f64 = theta * theta;
        let (v, chi_delta) = sim3_v(sigma, z);
        let chi = chi_delta / (sigma * sigma + z);
        let direct = v[0] - z * v[2];
        assert!(
            (chi - direct).abs() <= 8.0 * u * v[0].abs().max(z * v[2]),
            "({sigma}, {theta}): {chi} against {direct}"
        );
    }
}
