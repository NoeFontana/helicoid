//! Aggregation per `(fn, stratum, precision, subject)` and the two renderings of the rows: the
//! result CSV (`docs/PHASE1.md` §5) and the table `just conformance` prints.

use std::collections::BTreeMap;
use std::fmt::Write;

use helicoid_linalg::Precision;

use super::metric::Score;

/// The header of the result CSV, exactly.
pub(crate) const HEADER: &str =
    "fn,stratum,precision,subject,subject_version,n,max_u,p99_u,argmax_id,nonfinite,git_rev";

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Row {
    pub(crate) fn_id: String,
    pub(crate) stratum: String,
    pub(crate) precision: Precision,
    pub(crate) subject: String,
    pub(crate) version: String,
    /// Records evaluated, non-finite and unscored ones included.
    pub(crate) n: usize,
    /// The largest error of the scored records; NaN when there are none (`NaN` in the CSV).
    pub(crate) max_u: f64,
    /// Nearest-rank p99 of the same records: the `⌈0.99·m⌉`-th smallest of `m`. At `m = 64` it is
    /// the max; a stratum shows it apart from the max from `m = 100` on. NaN as `max_u` is.
    pub(crate) p99_u: f64,
    /// The first record with `max_u`; with none scored, the first non-finite record, else the
    /// stratum's first record.
    pub(crate) argmax_id: u64,
    pub(crate) nonfinite: usize,
}

/// One function's scores of one subject, grouped by stratum in order of first appearance.
#[derive(Default)]
pub(crate) struct Aggregate {
    groups: Vec<Group>,
    index: BTreeMap<String, usize>,
}

#[derive(Default)]
struct Group {
    stratum: String,
    n: usize,
    first_id: u64,
    finite: Vec<(u64, f64)>,
    nonfinite: usize,
    first_nonfinite: u64,
}

impl Aggregate {
    pub(crate) fn add(&mut self, id: u64, stratum: &str, score: Score) {
        let at = *self.index.entry(stratum.to_string()).or_insert_with(|| {
            self.groups.push(Group {
                stratum: stratum.to_string(),
                first_id: id,
                ..Group::default()
            });
            self.groups.len() - 1
        });
        let g = &mut self.groups[at];
        g.n += 1;
        match score {
            Score::Finite(f) => g.finite.push((id, f)),
            Score::NonFinite => {
                if g.nonfinite == 0 {
                    g.first_nonfinite = id;
                }
                g.nonfinite += 1;
            }
            Score::Unscored => {}
        }
    }

    /// One row per stratum, in order of first appearance.
    pub(crate) fn rows(
        self,
        fn_id: &str,
        precision: Precision,
        subject: &str,
        version: &str,
    ) -> Vec<Row> {
        self.groups
            .into_iter()
            .map(|g| {
                let mut sorted: Vec<f64> = g.finite.iter().map(|&(_, f)| f).collect();
                sorted.sort_by(f64::total_cmp);
                let p99_u = match sorted.len() {
                    0 => f64::NAN,
                    m => sorted[(99 * m).div_ceil(100) - 1],
                };
                let worst = g
                    .finite
                    .iter()
                    .fold(None, |best: Option<(u64, f64)>, &(id, f)| {
                        best.filter(|&(_, m)| f <= m).or(Some((id, f)))
                    });
                Row {
                    fn_id: fn_id.to_string(),
                    stratum: g.stratum,
                    precision,
                    subject: subject.to_string(),
                    version: version.to_string(),
                    n: g.n,
                    max_u: worst.map_or(f64::NAN, |(_, f)| f),
                    p99_u,
                    argmax_id: worst.map_or(
                        if g.nonfinite > 0 {
                            g.first_nonfinite
                        } else {
                            g.first_id
                        },
                        |(id, _)| id,
                    ),
                    nonfinite: g.nonfinite,
                }
            })
            .collect()
    }
}

pub(crate) fn precision_name(p: Precision) -> &'static str {
    match p {
        Precision::F32 => "f32",
        Precision::F64 => "f64",
    }
}

/// The shortest decimal that reads back as the same binary64: identical bytes on every machine,
/// and no loss for the exact no-regress bar.
fn float(x: f64) -> String {
    format!("{x:e}")
}

fn quoted(s: &str) -> String {
    if s.contains([',', '"', '\n']) {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

/// The result CSV: [`HEADER`], then `rows` in the order given. Every score is finite or the cell
/// is `NaN` (a stratum with none scored); `inf` never appears, since an overflowing error is
/// non-finite (`metric::Score`).
pub(crate) fn csv(rows: &[Row], git_rev: &str) -> String {
    let mut out = format!("{HEADER}\n");
    for r in rows {
        let cells = [
            quoted(&r.fn_id),
            quoted(&r.stratum),
            precision_name(r.precision).to_string(),
            quoted(&r.subject),
            quoted(&r.version),
            r.n.to_string(),
            float(r.max_u),
            float(r.p99_u),
            r.argmax_id.to_string(),
            r.nonfinite.to_string(),
            quoted(git_rev),
        ];
        out.push_str(&cells.join(","));
        out.push('\n');
    }
    out
}

/// `rows` by `max_u` descending, ties in the order given, as aligned columns. A stratum with a
/// non-finite record and none scored comes first, one with nothing scored and nothing wrong last.
pub(crate) fn table(rows: &[Row]) -> String {
    let key = |r: &Row| {
        if r.max_u.is_nan() && r.nonfinite == 0 {
            f64::NEG_INFINITY
        } else {
            r.max_u
        }
    };
    let mut sorted: Vec<&Row> = rows.iter().collect();
    sorted.sort_by(|a, b| key(b).total_cmp(&key(a)));
    let header = [
        "fn",
        "stratum",
        "prec",
        "subject",
        "n",
        "max_u",
        "p99_u",
        "argmax",
        "nonfinite",
    ];
    let cells: Vec<[String; 9]> = sorted
        .iter()
        .map(|r| {
            [
                r.fn_id.clone(),
                r.stratum.clone(),
                precision_name(r.precision).to_string(),
                r.subject.clone(),
                r.n.to_string(),
                format!("{:.4e}", r.max_u),
                format!("{:.4e}", r.p99_u),
                r.argmax_id.to_string(),
                r.nonfinite.to_string(),
            ]
        })
        .collect();
    let width: Vec<usize> = (0..9)
        .map(|c| {
            cells
                .iter()
                .map(|r| r[c].len())
                .fold(header[c].len(), usize::max)
        })
        .collect();
    let mut out = String::new();
    let mut line = |cols: &mut dyn Iterator<Item = &str>| {
        for (c, text) in cols.enumerate() {
            let w = width[c];
            let _ = if c < 4 {
                write!(out, "{text:<w$}  ")
            } else {
                write!(out, "{text:>w$}  ")
            };
        }
        out.truncate(out.trim_end().len());
        out.push('\n');
    };
    line(&mut header.into_iter());
    for r in &cells {
        line(&mut r.iter().map(String::as_str));
    }
    out
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    fn agg(scores: &[(&str, Score)]) -> Vec<Row> {
        let mut a = Aggregate::default();
        for (i, (stratum, s)) in scores.iter().enumerate() {
            a.add(i as u64, stratum, *s);
        }
        a.rows("f", Precision::F64, "s", "v")
    }

    #[test]
    fn nearest_rank_p99_is_the_max_below_100_records() {
        let scores = |m: usize| -> Vec<(&'static str, Score)> {
            (0..m)
                .map(|i| ("s", Score::Finite(((i * 37) % m) as f64)))
                .collect()
        };
        let p99 = |m: usize| agg(&scores(m))[0].p99_u;
        assert_eq!(p99(1), 0.0);
        assert_eq!(p99(64), 63.0);
        assert_eq!(p99(99), 98.0);
        // ⌈0.99·100⌉ = 99: the second largest; ⌈0.99·101⌉ = 100; ⌈0.99·801⌉ = 793.
        assert_eq!(p99(100), 98.0);
        assert_eq!(p99(101), 99.0);
        assert_eq!(p99(801), 792.0);
    }

    #[test]
    fn strata_keep_first_appearance_order_and_ties_go_to_the_first_record() {
        let rows = agg(&[
            ("b", Score::Finite(2.0)),
            ("a", Score::Finite(1.0)),
            ("b", Score::Finite(5.0)),
            ("a", Score::Finite(1.0)),
            ("b", Score::Finite(5.0)),
        ]);
        let got: Vec<_> = rows
            .iter()
            .map(|r| (r.stratum.as_str(), r.n, r.max_u, r.argmax_id))
            .collect();
        assert_eq!(got, [("b", 3, 5.0, 2), ("a", 2, 1.0, 1)]);
    }

    #[test]
    fn non_finite_records_are_counted_apart_from_the_max_and_the_p99() {
        let rows = agg(&[
            ("s", Score::Finite(3.0)),
            ("s", Score::NonFinite),
            ("s", Score::Finite(1.0)),
            ("t", Score::NonFinite),
            ("t", Score::NonFinite),
        ]);
        assert_eq!(
            (rows[0].n, rows[0].nonfinite, rows[0].max_u, rows[0].p99_u),
            (3, 1, 3.0, 3.0)
        );
        assert_eq!((rows[1].n, rows[1].nonfinite, rows[1].argmax_id), (2, 2, 3));
        assert!(rows[1].max_u.is_nan() && rows[1].p99_u.is_nan());
    }

    #[test]
    fn unscored_records_are_counted_and_name_their_first_record() {
        let rows = agg(&[
            ("s", Score::Unscored),
            ("s", Score::Unscored),
            ("t", Score::Unscored),
            ("t", Score::NonFinite),
            ("u", Score::Unscored),
            ("u", Score::Finite(4.0)),
        ]);
        let got: Vec<_> = rows
            .iter()
            .map(|r| (r.n, r.nonfinite, r.max_u.is_nan(), r.argmax_id))
            .collect();
        assert_eq!(got, [(2, 0, true, 0), (2, 1, true, 3), (2, 0, false, 5)]);
    }

    #[test]
    fn csv_has_the_exact_header_and_shortest_round_trip_floats() {
        let rows = agg(&[
            ("theta:1e-3", Score::Finite(0.5)),
            ("theta:1e-3", Score::Finite(1.0 / 3.0)),
        ]);
        let text = csv(&rows, "abc123");
        let want = "fn,stratum,precision,subject,subject_version,n,max_u,p99_u,argmax_id,nonfinite,git_rev\n\
                    f,theta:1e-3,f64,s,v,2,5e-1,5e-1,0,0,abc123\n";
        assert_eq!(text, want);
        let nan = agg(&[("s", Score::NonFinite)]);
        assert!(csv(&nan, "r").contains(",NaN,NaN,0,1,r\n"));
        let odd = Row {
            stratum: "a,\"b\"".into(),
            ..rows[0].clone()
        };
        assert!(csv(&[odd], "r").contains("\"a,\"\"b\"\"\""));
        assert_eq!(float(1.0 / 3.0), "3.333333333333333e-1");
    }

    #[test]
    fn the_table_is_sorted_by_max_descending_and_stable() {
        let rows = agg(&[
            ("a", Score::Finite(1.0)),
            ("b", Score::Finite(9.0)),
            ("c", Score::Finite(1.0)),
            ("d", Score::NonFinite),
        ]);
        let text = table(&rows);
        let order: Vec<&str> = text
            .lines()
            .skip(1)
            .map(|l| l.split_whitespace().nth(1).unwrap_or(""))
            .collect();
        assert_eq!(order, ["d", "b", "a", "c"]);
        // Nothing scored and nothing wrong sorts last, under every finite score.
        let rows = agg(&[("a", Score::Unscored), ("b", Score::Finite(0.0))]);
        assert!(table(&rows)
            .lines()
            .nth(1)
            .is_some_and(|l| l.contains(" b ")));
        assert!(text
            .lines()
            .next()
            .is_some_and(|h| h.starts_with("fn ") && h.ends_with("nonfinite")));
        assert!(text.contains("9.0000e0"));
    }
}
