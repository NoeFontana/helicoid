//! The harness against its own oracles, over the committed corpus: a subject that returns the
//! correctly rounded reference, and one that returns the neighbouring ulp. Both know the answer,
//! so a wrong metric, floor or reader shows here before any real subject exists.

use helicoid_linalg::Precision;

use super::corpus::Record;
use super::metric::{self, Rule};
use super::report::Row;
use super::subject::{Output, Registered, Subject};
use super::testkit::{Fixed, Perfect};
use super::{corpus, corpus_dir, evaluate};

/// The rows of `subject` over the corpus files `keep` selects.
fn rows(subject: impl Subject + 'static, keep: impl Fn(&str) -> bool) -> Result<Vec<Row>, String> {
    let dir = corpus_dir()?;
    let entries: Vec<_> = corpus::manifest(&dir)?
        .into_iter()
        .filter(|e| keep(&e.fn_id))
        .collect();
    let registered = [Registered {
        version: "test".into(),
        subject: Box::new(subject),
    }];
    Ok(evaluate(&dir, &entries, &registered, Precision::F64)?.remove(0))
}

/// Ids whose every output field is one number: the ones where the score of a rounding is bounded
/// by `u` per record whatever the floor (a vector of subnormals can read up to `√n`).
fn scalar_output(fn_id: &str) -> Result<bool, String> {
    let dir = corpus_dir()?;
    let text =
        std::fs::read_to_string(dir.join(format!("{fn_id}.jsonl"))).map_err(|e| e.to_string())?;
    let first = corpus::parse_line(text.lines().next().unwrap_or_default())?;
    Ok(first.reference.values().all(|t| t.data.len() == 1))
}

fn scalar_ids() -> Result<Vec<String>, String> {
    let mut ids = Vec::new();
    for e in corpus::manifest(&corpus_dir()?)? {
        if scalar_output(&e.fn_id)? {
            ids.push(e.fn_id);
        }
    }
    Ok(ids)
}

#[test]
fn every_family_has_a_rule_and_the_corpus_has_exactly_the_fields_it_scores() -> Result<(), String> {
    let dir = corpus_dir()?;
    for e in corpus::manifest(&dir)? {
        let rule = metric::rule(&e.fn_id).ok_or(format!("{} has no rule", e.fn_id))?;
        let text = std::fs::read_to_string(dir.join(format!("{}.jsonl", e.fn_id)))
            .map_err(|e| e.to_string())?;
        let first = corpus::parse_line(text.lines().next().unwrap_or_default())?;
        if let Rule::Forward(fields) = rule {
            let mut scored: Vec<&str> = fields.iter().map(|f| f.field).collect();
            scored.sort_unstable();
            let corpus: Vec<&str> = first.reference.keys().map(String::as_str).collect();
            assert_eq!(scored, corpus, "{}", e.fn_id);
        }
    }
    Ok(())
}

#[test]
fn a_perfectly_rounded_subject_scores_at_most_one_on_every_scalar_id_and_stratum(
) -> Result<(), String> {
    let ids = scalar_ids()?;
    assert!(
        ids.iter().any(|i| i == "coeff_e") && ids.iter().any(|i| i == "so2_log"),
        "{ids:?}"
    );
    let rows = rows(Perfect::exact(), |id| ids.iter().any(|i| i == id))?;
    assert!(rows.len() > 200, "{} rows", rows.len());
    for r in &rows {
        assert!(
            r.max_u <= 1.0 && r.nonfinite == 0,
            "{} {}: {}",
            r.fn_id,
            r.stratum,
            r.max_u
        );
    }
    Ok(())
}

#[test]
fn the_neighbouring_ulp_scores_above_one_half_and_at_most_three_on_every_scalar_id_and_stratum(
) -> Result<(), String> {
    let ids = scalar_ids()?;
    let rows = rows(Perfect::next_up(), |id| ids.iter().any(|i| i == id))?;
    for r in &rows {
        assert!(
            r.max_u > 0.5 && r.max_u <= 3.0,
            "{} {}: {}",
            r.fn_id,
            r.stratum,
            r.max_u
        );
    }
    Ok(())
}

#[test]
fn a_perfectly_rounded_subject_is_finite_and_within_two_u_on_every_forward_id() -> Result<(), String>
{
    let rows = rows(Perfect::exact(), |_| true)?;
    assert!(rows.len() > 1900, "{} rows", rows.len());
    for r in &rows {
        let backward_only = r.fn_id == "so3_from_matrix";
        assert!(
            r.nonfinite == 0 && (backward_only || r.max_u <= 2.0),
            "{} {}: {}",
            r.fn_id,
            r.stratum,
            r.max_u
        );
    }
    // Every record reaches an aggregate, and a backward-only id scores nothing but is counted.
    for e in corpus::manifest(&corpus_dir()?)? {
        let mine = rows.iter().filter(|r| r.fn_id == e.fn_id);
        assert_eq!(mine.map(|r| r.n).sum::<usize>(), e.records, "{}", e.fn_id);
    }
    let unscored = |r: &&Row| r.fn_id == "so3_from_matrix";
    assert!(rows.iter().filter(unscored).count() > 5);
    assert!(rows
        .iter()
        .filter(unscored)
        .all(|r| r.max_u.is_nan() && r.p99_u.is_nan()));
    Ok(())
}

/// The correctly rounded reference with every value of `field` negated.
fn negated(field: &'static str) -> impl Fn(&Record) -> Output {
    move |rec| {
        let mut out = Perfect::exact().eval("", rec, Precision::F64);
        if let Some(values) = out.get_mut(field) {
            values.iter_mut().for_each(|x| *x = -*x);
        }
        out
    }
}

#[test]
fn so3_log_forgives_a_sign_only_at_w_plus_zero_and_sen3_log_never() -> Result<(), String> {
    let so3 = rows(Fixed::new("flip", negated("phi")), |id| id == "so3_log")?;
    for r in &so3 {
        match r.stratum.as_str() {
            "q:w0" => assert!(r.max_u <= 1.0, "{}", r.max_u),
            "theta:exact0" => assert!(r.max_u <= 1.0, "{}", r.max_u),
            s => assert!(r.max_u > 1e14, "{s}: {}", r.max_u),
        }
    }
    assert!(so3.iter().any(|r| r.stratum == "q:w0"));
    let sen3 = rows(Fixed::new("flip", negated("tau")), |id| id == "sen3_log_n1")?;
    let w0 = sen3.iter().find(|r| r.stratum == "q:w0").ok_or("q:w0")?;
    assert!(w0.max_u > 1e15, "{}", w0.max_u);
    Ok(())
}

#[test]
fn a_gross_error_against_a_zero_reference_is_non_finite_not_infinite() -> Result<(), String> {
    let gross = Fixed::new("gross", |rec| {
        let mut out = Perfect::exact().eval("", rec, Precision::F64);
        out.insert("phi".into(), vec![1.0, 0.0, 0.0]);
        out
    });
    let rows = rows(gross, |id| id == "so3_log")?;
    let zero = rows
        .iter()
        .find(|r| r.stratum == "theta:exact0")
        .ok_or("theta:exact0")?;
    assert!(zero.n > 0 && zero.nonfinite == zero.n, "{zero:?}");
    assert!(rows.iter().all(|r| !r.max_u.is_infinite()));
    Ok(())
}
