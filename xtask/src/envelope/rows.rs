//! Reading what the bars compare: a subject's result CSV (`conformance::report::csv` writes it),
//! which is also the format of the committed baseline, and the key rows are paired by.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use helicoid_linalg::Precision;

use crate::conformance::report::{precision_name, Row, HEADER};

use super::bars::is_scored;

/// `(fn, stratum, precision)`: what a bar is per (`docs/PHASE1.md` §8).
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(super) struct Key {
    fn_id: String,
    stratum: String,
    precision: &'static str,
}

impl Key {
    pub(super) fn of(r: &Row) -> Self {
        Self {
            fn_id: r.fn_id.clone(),
            stratum: r.stratum.clone(),
            precision: precision_name(r.precision),
        }
    }
}

impl fmt::Display for Key {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}/{}/{}", self.fn_id, self.stratum, self.precision)
    }
}

/// One subject's rows, in file order, of one name and one version.
pub(super) struct Subject {
    pub(super) name: String,
    pub(super) version: String,
    pub(super) rows: Vec<Row>,
}

/// What a baseline holds per key: the records the maximum is over, and the maximum.
pub(super) type Baseline = BTreeMap<Key, (usize, f64)>;

/// One CSV line's cells; a cell in quotes holds `,` and `""` for `"`, as `report::quoted` writes.
fn cells(line: &str) -> Result<Vec<String>, String> {
    let (mut out, mut cell, mut quoted) = (Vec::new(), String::new(), false);
    let mut chars = line.chars().peekable();
    while let Some(c) = chars.next() {
        match (c, quoted) {
            ('"', true) if chars.peek() == Some(&'"') => {
                cell.push('"');
                chars.next();
            }
            ('"', true) => quoted = false,
            ('"', false) if cell.is_empty() => quoted = true,
            (',', false) => out.push(std::mem::take(&mut cell)),
            (c, _) => cell.push(c),
        }
    }
    match quoted {
        true => Err("a quote is not closed".into()),
        false => {
            out.push(cell);
            Ok(out)
        }
    }
}

fn count<T: std::str::FromStr>(what: &str, s: &str) -> Result<T, String> {
    s.parse()
        .map_err(|_| format!("{what} `{s}` is not a count"))
}

/// An error in units of `u`, or `NaN` for a stratum nothing was scored in; never `inf` or negative.
fn error(what: &str, s: &str) -> Result<f64, String> {
    let x: f64 = s
        .parse()
        .map_err(|_| format!("{what} `{s}` is not a number"))?;
    match x.is_nan() || (x.is_finite() && x >= 0.0) {
        true => Ok(x),
        false => Err(format!("{what} `{s}` is neither NaN nor a finite error")),
    }
}

fn precision(s: &str) -> Result<Precision, String> {
    match s {
        "f64" => Ok(Precision::F64),
        "f32" => Ok(Precision::F32),
        _ => Err(format!("precision `{s}` is neither f64 nor f32")),
    }
}

fn row(line: &str) -> Result<Row, String> {
    let all = cells(line)?;
    let n = all.len();
    let Ok([fn_id, stratum, prec, subject, version, records, max_u, p99_u, argmax, nonfinite, _]) =
        <[String; 11]>::try_from(all)
    else {
        return Err(format!("{n} cells, expected 11"));
    };
    Ok(Row {
        precision: precision(&prec)?,
        n: count("n", &records)?,
        max_u: error("max_u", &max_u)?,
        p99_u: error("p99_u", &p99_u)?,
        argmax_id: count("argmax_id", &argmax)?,
        nonfinite: count("nonfinite", &nonfinite)?,
        fn_id,
        stratum,
        subject,
        version,
    })
}

/// The rows of `name`'s result CSV. A repeated key, another subject's row, two versions in one
/// file and a header that is not `report::HEADER` are errors: each would make a bar read a run
/// that is not one run of one subject.
pub(super) fn parse_results(name: &str, text: &str) -> Result<Subject, String> {
    let mut lines = text.lines();
    if lines.next() != Some(HEADER) {
        return Err(format!("the header is not `{HEADER}`"));
    }
    let (mut rows, mut seen) = (Vec::new(), BTreeSet::new());
    for (i, line) in lines.enumerate() {
        let at = |e: String| format!("line {}: {e}", i + 2);
        let r = row(line).map_err(at)?;
        if r.subject != name {
            return Err(at(format!("subject `{}`, expected `{name}`", r.subject)));
        }
        if !seen.insert(Key::of(&r)) {
            return Err(at(format!("`{}` twice", Key::of(&r))));
        }
        rows.push(r);
    }
    let version = rows.first().map(|r| r.version.clone()).unwrap_or_default();
    if let Some(r) = rows.iter().find(|r| r.version != version) {
        return Err(format!(
            "versions `{version}` and `{}` in one file",
            r.version
        ));
    }
    Ok(Subject {
        name: name.to_string(),
        version,
        rows,
    })
}

/// The baseline `s` makes: per key the records and the maximum of each row. A row with no maximum
/// (non-finite, or nothing scored) is no bar, so it is no baseline row.
pub(super) fn baseline_of(s: &Subject) -> Result<Baseline, String> {
    let mut out = Baseline::new();
    for r in &s.rows {
        if !is_scored(r) {
            return Err(format!("`{}` has no maximum", Key::of(r)));
        }
        out.insert(Key::of(r), (r.n, r.max_u));
    }
    Ok(out)
}

#[cfg(test)]
pub(super) mod testing {
    use super::*;

    /// A scored row of 64 records with `max_u` and `p99_u` half of it, so that a column that reads
    /// one for the other shows.
    pub(crate) fn row(subject: &str, fn_id: &str, stratum: &str, max_u: f64) -> Row {
        Row {
            fn_id: fn_id.into(),
            stratum: stratum.into(),
            precision: Precision::F64,
            subject: subject.into(),
            version: format!("{subject}@1"),
            n: 64,
            max_u,
            p99_u: max_u / 2.0,
            argmax_id: 0,
            nonfinite: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::testing::row;
    use super::*;
    use crate::conformance::report::csv;

    #[test]
    fn what_the_harness_writes_reads_back_exactly() -> Result<(), String> {
        let odd = row("s", "f", "a,\"b\"", 1.0 / 3.0);
        let nan = Row {
            max_u: f64::NAN,
            p99_u: f64::NAN,
            nonfinite: 2,
            ..row("s", "g", "z", 0.0)
        };
        let f32_row = Row {
            precision: Precision::F32,
            ..row("s", "f", "z", 5e-324)
        };
        let rows = vec![odd, row("s", "f", "z", 2.5), nan, f32_row];
        let read = parse_results("s", &csv(&rows, "rev, with comma"))?;
        assert_eq!((read.name.as_str(), read.version.as_str()), ("s", "s@1"));
        assert_eq!(read.rows.len(), 4);
        for (a, b) in rows.iter().zip(&read.rows) {
            let same = |x: f64, y: f64| x.to_bits() == y.to_bits();
            assert!(same(a.max_u, b.max_u) && same(a.p99_u, b.p99_u), "{a:?}");
            assert_eq!(
                (&a.stratum, a.precision, a.n),
                (&b.stratum, b.precision, b.n)
            );
            assert_eq!((a.nonfinite, a.argmax_id), (b.nonfinite, b.argmax_id));
        }
        Ok(())
    }

    #[test]
    fn a_file_that_is_not_one_run_of_one_subject_is_an_error() {
        let good = row("s", "f", "z", 1.0);
        let err = |name: &str, text: String| parse_results(name, &text).err().unwrap_or_default();
        let one = csv(std::slice::from_ref(&good), "r");
        assert!(err("t", one.clone()).contains("subject `s`, expected `t`"));
        assert!(err("s", one.replace("fn,", "fun,")).contains("header"));
        let twice = csv(&[good.clone(), good.clone()], "r");
        assert!(err("s", twice).contains("line 3: `f/z/f64` twice"));
        let other = Row {
            version: "s@2".into(),
            stratum: "y".into(),
            ..good.clone()
        };
        assert!(err("s", csv(&[good, other], "r")).contains("versions `s@1` and `s@2`"));
        for (from, to, needle) in [
            (",1e0,", ",inf,", "neither NaN nor a finite error"),
            (",1e0,", ",-1e0,", "neither NaN nor a finite error"),
            (",1e0,", ",x,", "not a number"),
            (",f64,", ",f16,", "neither f64 nor f32"),
            (",64,", ",-1,", "not a count"),
            ("z,f64", "\"z,f64", "not closed"),
        ] {
            let text = csv(&[row("s", "f", "z", 1.0)], "r").replacen(from, to, 1);
            assert!(err("s", text.clone()).contains(needle), "{text}");
        }
        assert!(err("s", format!("{HEADER}\na,b\n")).contains("2 cells, expected 11"));
    }

    #[test]
    fn a_baseline_is_the_scored_rows_of_a_subject() -> Result<(), String> {
        let rows = vec![row("s", "f", "a", 0.1), row("s", "f", "z", 5e-324)];
        let subject = Subject {
            name: "s".into(),
            version: "s@1".into(),
            rows,
        };
        let b = baseline_of(&subject)?;
        let keys: Vec<String> = b.keys().map(ToString::to_string).collect();
        assert_eq!(keys, ["f/a/f64", "f/z/f64"]);
        assert_eq!(b[&Key::of(&subject.rows[1])], (64, 5e-324));
        let mut bad = subject;
        bad.rows[0].max_u = f64::NAN;
        assert!(baseline_of(&bad)
            .err()
            .unwrap_or_default()
            .contains("no maximum"));
        Ok(())
    }
}
