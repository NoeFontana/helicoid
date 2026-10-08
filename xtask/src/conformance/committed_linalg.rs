//! `0056` decision 5: the committed maximum of every stratum of the routines D7 did not reach,
//! reproduced bit for bit by the shipped subject.
//!
//! No-regress without the envelope's first bless, and in both directions: an improvement is a
//! re-recorded row as a regression is a failure. CI runs it on x86_64 and aarch64, so each row is
//! also D16's bit identity across the two. `measure_linalg` prints the tables below; paste them
//! over [`F64`] and [`F32`] when a change to one of these routines is meant.

use helicoid_linalg::Precision;

use super::{corpus, corpus_dir, evaluate};
use crate::shipped;

/// The ids of `0056` decision 1.
const IDS: [&str; 13] = [
    "solve_cubic",
    "eig3",
    "chol_n3",
    "chol_n6",
    "chol_solve_n3",
    "chol_solve_n6",
    "quat_renormalize",
    "real_sqrt",
    "real_cbrt",
    "real_sin_cos",
    "real_acos",
    "real_atan2",
    "real_div",
];

/// `(fn, stratum, max_u)` of the binary64 strata.
const F64: &[(&str, &str, f64)] = &[
    ("chol_n3", "chol:spd", 0.8946325689466859),
    ("chol_n3", "chol:cond-1e-4", 55.58068627719461),
    ("chol_n3", "chol:cond-1e-8", 26060.404622832248),
    ("chol_n3", "chol:cond-1e-12", 1834397.862150613),
    ("chol_n3", "chol:diag-scale", 1.1079982581157728),
    ("chol_n3", "chol:indefinite", 0.0),
    ("chol_n6", "chol:spd", 0.7039317735237633),
    ("chol_n6", "chol:cond-1e-4", 61.01493814619754),
    ("chol_n6", "chol:cond-1e-8", 4249.786823495965),
    ("chol_n6", "chol:cond-1e-12", 2593879.6757169846),
    ("chol_n6", "chol:diag-scale", 1.1465717134477733),
    ("chol_n6", "chol:indefinite", 0.0),
    ("chol_solve_n3", "chol:spd", 2.8148682219084407),
    ("chol_solve_n3", "chol:cond-1e-4", 2680.7736852254807),
    ("chol_solve_n3", "chol:cond-1e-8", 63332260.698645875),
    ("chol_solve_n3", "chol:cond-1e-12", 303254211025.0341),
    ("chol_solve_n3", "chol:diag-scale", 1.8783310159551247),
    ("chol_solve_n6", "chol:spd", 2.7425717175142914),
    ("chol_solve_n6", "chol:cond-1e-4", 6017.081512624412),
    ("chol_solve_n6", "chol:cond-1e-8", 42515988.89538545),
    ("chol_solve_n6", "chol:cond-1e-12", 310101739113.7394),
    ("chol_solve_n6", "chol:diag-scale", 4.592010440093427),
    ("eig3", "eig:random", 22.69637423897279),
    ("eig3", "eig:gap-1e-0/bottom", 5.632318004600339),
    ("eig3", "eig:gap-1e-1/bottom", 4.5988606812250135),
    ("eig3", "eig:gap-1e-2/bottom", 47.92107028542843),
    ("eig3", "eig:gap-1e-3/bottom", 326.5960597488452),
    ("eig3", "eig:gap-1e-4/bottom", 4259.333913046709),
    ("eig3", "eig:gap-1e-5/bottom", 39449.30564745859),
    ("eig3", "eig:gap-1e-6/bottom", 537279.9877232165),
    ("eig3", "eig:gap-1e-7/bottom", 3707546.7665354274),
    ("eig3", "eig:gap-1e-8/bottom", 100913902.66343957),
    ("eig3", "eig:gap-1e-9/bottom", 53853339.43977939),
    ("eig3", "eig:gap-1e-10/bottom", 58822649.89179961),
    ("eig3", "eig:gap-1e-11/bottom", 71359868.47346932),
    ("eig3", "eig:gap-1e-12/bottom", 84220936.55792174),
    ("eig3", "eig:gap-1e-0/top", 5.475903782991654),
    ("eig3", "eig:gap-1e-1/top", 4.94627360742461),
    ("eig3", "eig:gap-1e-2/top", 55.770645933131945),
    ("eig3", "eig:gap-1e-3/top", 391.96738312022046),
    ("eig3", "eig:gap-1e-4/top", 2905.2403743816553),
    ("eig3", "eig:gap-1e-5/top", 29090.98255903702),
    ("eig3", "eig:gap-1e-6/top", 378775.69797502796),
    ("eig3", "eig:gap-1e-7/top", 5925471.788193611),
    ("eig3", "eig:gap-1e-8/top", 107122939.70572017),
    ("eig3", "eig:gap-1e-9/top", 68007307.42535928),
    ("eig3", "eig:gap-1e-10/top", 70164435.52491628),
    ("eig3", "eig:gap-1e-11/top", 52145618.66884013),
    ("eig3", "eig:gap-1e-12/top", 63197835.28545676),
    ("eig3", "eig:triple", 2.2554304411603985),
    ("eig3", "eig:rank1", 63270842.984606355),
    ("eig3", "eig:scale-up", 22.511124803937662),
    ("eig3", "eig:scale-down", 15.461887691971556),
    ("quat_renormalize", "renorm:eta-2^-27", 1.4113137658448696),
    ("quat_renormalize", "renorm:eta-2^-30", 1.562244977607657),
    ("quat_renormalize", "renorm:eta-2^-40", 1.2309318790786339),
    ("quat_renormalize", "renorm:eta-2^-52", 1.418197746297618),
    ("quat_renormalize", "renorm:eta-edge", 2.054067917991086),
    ("real_acos", "x:interior", 1.4292739672947972),
    ("real_acos", "x:near+1", 1.36750190045526),
    ("real_acos", "x:near-1", 1.36750190045526),
    ("real_acos", "x:tiny", 0.3511161042472099),
    ("real_atan2", "yx:generic", 2.2048315342842244),
    ("real_atan2", "yx:ratio-1e-8", 2.0468900050436),
    ("real_atan2", "yx:ratio-1e-100", 1.5844965782660128),
    ("real_atan2", "yx:ratio-1e-300", 1.5520307300777068),
    ("real_atan2", "yx:ratio-1e+8", 1.9660061781056537),
    ("real_atan2", "yx:ratio-1e+100", 1.5130407237361612),
    ("real_atan2", "yx:ratio-1e+300", 1.4519165423436347),
    ("real_cbrt", "x:1e-300", 2.7449087858008165),
    ("real_cbrt", "x:1e-100", 2.431719345997256),
    ("real_cbrt", "x:1e-10", 2.237978111023069),
    ("real_cbrt", "x:1e-1", 3.4658482887582416),
    ("real_cbrt", "x:1e0", 3.2217936233858016),
    ("real_cbrt", "x:1e1", 3.1550994879215755),
    ("real_cbrt", "x:1e10", 2.90242394611477),
    ("real_cbrt", "x:1e100", 3.144910884910384),
    ("real_cbrt", "x:1e300", 2.2214861294543504),
    ("real_cbrt", "x:subnormal", 2.6404449005484603),
    ("real_div", "nd:generic", 1.4167055841482492),
    ("real_div", "nd:wide", 1.205283488858324),
    ("real_sin_cos", "x:tiny", 8.624328434579866e-8),
    ("real_sin_cos", "x:small", 0.8560314598447317),
    ("real_sin_cos", "x:moderate", 1.019543605711025),
    ("real_sin_cos", "x:large", 0.9242752035959709),
    ("real_sin_cos", "x:near-k-pi/2", 0.8369921005849935),
    ("real_sqrt", "x:1e-300", 1.4093799217226632),
    ("real_sqrt", "x:1e-100", 1.2915590393964518),
    ("real_sqrt", "x:1e-10", 1.3545740271299138),
    ("real_sqrt", "x:1e-1", 1.234018402081571),
    ("real_sqrt", "x:1e0", 1.1822493215506624),
    ("real_sqrt", "x:1e1", 1.16547184528416),
    ("real_sqrt", "x:1e10", 1.1680344221155263),
    ("real_sqrt", "x:1e100", 1.341890011009563),
    ("real_sqrt", "x:1e300", 1.3346137626205352),
    ("real_sqrt", "x:subnormal", 1.2752358462221747),
    ("solve_cubic", "cubic:distinct", 37.0561032994736),
    ("solve_cubic", "cubic:double", 154457280.5870304),
    ("solve_cubic", "cubic:triple", 0.0),
    ("solve_cubic", "cubic:one-real", 1.7001805707579358),
    ("solve_cubic", "cubic:near-double-1e-2", 509.27880489346023),
    ("solve_cubic", "cubic:near-double-1e-4", 108300.5935573035),
    ("solve_cubic", "cubic:near-double-1e-6", 15878262.95250792),
    ("solve_cubic", "cubic:near-double-1e-8", 90869871.39177251),
    ("solve_cubic", "cubic:one-real-p-small", 1.7102572974660597),
    ("solve_cubic", "cubic:coeff-scale-up", 16.45779949292739),
    ("solve_cubic", "cubic:coeff-scale-down", 51.14647623707047),
    (
        "solve_cubic",
        "cubic:one-real-small-root",
        0.0054587359395373595,
    ),
];

/// `(fn, stratum, max_u)` of the `@f32` strata.
const F32: &[(&str, &str, f64)] = &[
    ("chol_n3", "chol:spd@f32", 0.851502723047995),
    ("chol_n3", "chol:cond-1e-2@f32", 4.964220695762808),
    ("chol_n3", "chol:cond-1e-4@f32", 112.49970195975011),
    ("chol_n3", "chol:cond-1e-6@f32", 1038.5534527243672),
    ("chol_n3", "chol:diag-scale@f32", 1.066930345401299),
    ("chol_n3", "chol:indefinite@f32", 0.0),
    ("chol_n6", "chol:spd@f32", 0.7975137437323015),
    ("chol_n6", "chol:cond-1e-2@f32", 3.999870973941397),
    ("chol_n6", "chol:cond-1e-4@f32", 62.188412048040696),
    ("chol_n6", "chol:cond-1e-6@f32", 474.10244743837626),
    ("chol_n6", "chol:diag-scale@f32", 1.0940822496240288),
    ("chol_n6", "chol:indefinite@f32", 0.0),
    ("chol_solve_n3", "chol:spd@f32", 1.933808030458939),
    ("chol_solve_n3", "chol:cond-1e-2@f32", 54.17610456448606),
    ("chol_solve_n3", "chol:cond-1e-4@f32", 5152.468014823599),
    ("chol_solve_n3", "chol:cond-1e-6@f32", 504739.8109168917),
    ("chol_solve_n3", "chol:diag-scale@f32", 3.1909776235293412),
    ("chol_solve_n6", "chol:spd@f32", 3.503763074634495),
    ("chol_solve_n6", "chol:cond-1e-2@f32", 31.620244174736467),
    ("chol_solve_n6", "chol:cond-1e-4@f32", 4783.657281312161),
    ("chol_solve_n6", "chol:cond-1e-6@f32", 369498.1618248905),
    ("chol_solve_n6", "chol:diag-scale@f32", 3.064501676959737),
    ("eig3", "eig:random@f32", 6.273804654657301),
    ("eig3", "eig:gap-1e-0/bottom@f32", 5.245015331866646),
    ("eig3", "eig:gap-1e-1/bottom@f32", 5.937185292228242),
    ("eig3", "eig:gap-1e-2/bottom@f32", 55.275513769397364),
    ("eig3", "eig:gap-1e-3/bottom@f32", 243.58853981717897),
    ("eig3", "eig:gap-1e-4/bottom@f32", 2268.8398050464043),
    ("eig3", "eig:gap-1e-5/bottom@f32", 2120.894403753257),
    ("eig3", "eig:gap-1e-6/bottom@f32", 2330.3722612227957),
    ("eig3", "eig:gap-1e-0/top@f32", 6.689091431635208),
    ("eig3", "eig:gap-1e-1/top@f32", 5.303287979069653),
    ("eig3", "eig:gap-1e-2/top@f32", 40.74883149820234),
    ("eig3", "eig:gap-1e-3/top@f32", 408.25876539677387),
    ("eig3", "eig:gap-1e-4/top@f32", 3682.8707709848195),
    ("eig3", "eig:gap-1e-5/top@f32", 2254.501731159185),
    ("eig3", "eig:gap-1e-6/top@f32", 3133.859477577778),
    ("eig3", "eig:triple@f32", 2.3793294302645696),
    ("eig3", "eig:rank1@f32", 3524.0701023503825),
    ("eig3", "eig:scale-up@f32", 11.950665709174782),
    ("eig3", "eig:scale-down@f32", 51.05561680873616),
    (
        "quat_renormalize",
        "renorm:eta-2^-12@f32",
        1.4716527222860658,
    ),
    (
        "quat_renormalize",
        "renorm:eta-2^-16@f32",
        1.1084966808592092,
    ),
    (
        "quat_renormalize",
        "renorm:eta-2^-20@f32",
        0.9747254412628126,
    ),
    (
        "quat_renormalize",
        "renorm:eta-2^-23@f32",
        1.3589338863391354,
    ),
    (
        "quat_renormalize",
        "renorm:eta-edge@f32",
        1.7527606454698232,
    ),
    ("real_acos", "x:interior@f32", 1.569892610758),
    ("real_acos", "x:near+1@f32", 1.4999999329447735),
    ("real_acos", "x:near-1@f32", 1.4999999329447735),
    ("real_acos", "x:tiny@f32", 1.9999999739192011),
    ("real_atan2", "yx:generic@f32", 1.9471729032796006),
    ("real_atan2", "yx:ratio-1e-4@f32", 1.4774569104626287),
    ("real_atan2", "yx:ratio-1e-15@f32", 1.5624183352687933),
    ("real_atan2", "yx:ratio-1e-30@f32", 1.3881208896636963),
    ("real_atan2", "yx:ratio-1e+4@f32", 1.453553547253106),
    ("real_atan2", "yx:ratio-1e+15@f32", 1.583339334476575),
    ("real_atan2", "yx:ratio-1e+30@f32", 1.3671131134033203),
    ("real_cbrt", "x:1e-37@f32", 3.6528471816191144),
    ("real_cbrt", "x:1e-10@f32", 2.5567921639279163),
    ("real_cbrt", "x:1e-1@f32", 2.4722057352498568),
    ("real_cbrt", "x:1e0@f32", 3.4097903968172565),
    ("real_cbrt", "x:1e1@f32", 2.521687995857178),
    ("real_cbrt", "x:1e10@f32", 3.2311844250306363),
    ("real_cbrt", "x:1e37@f32", 2.820464993199371),
    ("real_cbrt", "x:subnormal@f32", 2.4466267805138995),
    ("real_div", "nd:generic@f32", 1.329488952532492),
    ("real_div", "nd:wide@f32", 1.4416723978949668),
    ("real_sin_cos", "x:tiny@f32", 0.010809019401828107),
    ("real_sin_cos", "x:small@f32", 0.8967256490819028),
    ("real_sin_cos", "x:moderate@f32", 0.8976687850517485),
    ("real_sin_cos", "x:large@f32", 0.9218181361328055),
    ("real_sin_cos", "x:near-k-pi/2@f32", 0.9283828768601604),
    ("real_sqrt", "x:1e-37@f32", 1.2758175586676654),
    ("real_sqrt", "x:1e-10@f32", 1.279554236758641),
    ("real_sqrt", "x:1e-1@f32", 1.357725527045076),
    ("real_sqrt", "x:1e0@f32", 1.2926064541883875),
    ("real_sqrt", "x:1e1@f32", 1.278757052081233),
    ("real_sqrt", "x:1e10@f32", 1.398619799782542),
    ("real_sqrt", "x:1e37@f32", 1.3674029719326242),
    ("real_sqrt", "x:subnormal@f32", 1.0964127862704565),
    ("solve_cubic", "cubic:distinct@f32", 44.746992663515286),
    ("solve_cubic", "cubic:double@f32", 1809172.1779990862),
    ("solve_cubic", "cubic:triple@f32", 68493.4326087546),
    ("solve_cubic", "cubic:one-real@f32", 2.4550904125803092),
    (
        "solve_cubic",
        "cubic:near-double-1e-2@f32",
        1261.9892603896478,
    ),
    (
        "solve_cubic",
        "cubic:near-double-1e-4@f32",
        5666.5353639471505,
    ),
    (
        "solve_cubic",
        "cubic:near-double-1e-6@f32",
        5297.936203177212,
    ),
    (
        "solve_cubic",
        "cubic:near-double-1e-8@f32",
        5231.775054812178,
    ),
    (
        "solve_cubic",
        "cubic:one-real-p-small@f32",
        1.3888154836145314,
    ),
    ("solve_cubic", "cubic:coeff-scale-up@f32", 42.45508436606197),
    (
        "solve_cubic",
        "cubic:coeff-scale-down@f32",
        77.86191104182657,
    ),
    (
        "solve_cubic",
        "cubic:one-real-small-root@f32",
        0.06909626748909296,
    ),
];

/// The shipped subject's `(fn, stratum, max_u)` over [`IDS`] at `precision`.
fn measured(precision: Precision) -> Result<Vec<(String, String, f64)>, String> {
    let dir = corpus_dir()?;
    let entries: Vec<_> = corpus::manifest(&dir)?
        .into_iter()
        .filter(|e| IDS.contains(&e.fn_id.as_str()))
        .collect();
    let rows = evaluate(&dir, &entries, &[shipped::registered()], precision)?.remove(0);
    Ok(rows
        .into_iter()
        .map(|r| (r.fn_id, r.stratum, r.max_u))
        .collect())
}

/// `rows` as the Rust a table above is written in.
fn table(rows: &[(String, String, f64)]) -> String {
    rows.iter()
        .map(|(f, s, u)| format!("    ({f:?}, {s:?}, {u:?}),\n"))
        .collect()
}

#[test]
fn the_unverified_routines_reproduce_their_committed_rows() -> Result<(), String> {
    for (precision, committed) in [(Precision::F64, F64), (Precision::F32, F32)] {
        let got = measured(precision)?;
        let differs = |((f, s, u), (cf, cs, cu)): &(&(String, String, f64), &(&str, &str, f64))| {
            f != cf || s != cs || u.to_bits() != cu.to_bits()
        };
        let first = got.iter().zip(committed).find(differs);
        if got.len() != committed.len() || first.is_some() {
            return Err(format!(
                "{precision:?}: {} rows against {} committed, first difference {first:?}; \
                 if the change is meant, the measured table is:\n{}",
                got.len(),
                committed.len(),
                table(&got)
            ));
        }
    }
    Ok(())
}

#[test]
#[ignore = "a measurement, printed for 0056 and for the tables above"]
#[allow(clippy::print_stdout)]
fn measure_linalg() -> Result<(), String> {
    for precision in [Precision::F64, Precision::F32] {
        println!("{precision:?}:\n{}", table(&measured(precision)?));
    }
    Ok(())
}
