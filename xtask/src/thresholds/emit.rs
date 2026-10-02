//! `xtask/src/seeded/generated.rs` (`docs/PHASE1.md` §6, `0004` item 3): the seeded twin of
//! Phase 3's `coeffs/generated.rs`, a function of the sweep CSV and `coeff_series.jsonl` and of
//! nothing else, so `thresholds --check` compares it byte for byte. One `Switch<f64, terms>` per
//! coefficient, `K_F64` to `E_F64`, defined in `seeded::switch`.
//!
//! The file is compiled into this tool (`seeded:correct` runs it), so a hand edit that no longer
//! compiles stops `cargo xtask` from building, and neither `thresholds` nor `thresholds --check`
//! can name it: `git restore` the file first. One that compiles is named, file and line.
//!
//! Readings where §6 is silent, each the smallest (the first two are 0014 (draft) questions 12 and 13):
//!
//! - **The source revision** is two SHA-256 digests, of the CSV's bytes and of the series file's,
//!   which supplies the literals. A git revision cannot name the commit that contains the file and
//!   differs on every commit, which no byte-for-byte check survives.
//! - **The objective line** is one per constant, above it, from the CSV's own columns: §6's
//!   example shows one line for one constant, and a file holds six.
//! - **`terms`** is the array length, so the kernel's `terms` is `series.len()` and the constant
//!   holds the terms the sweep measured and no others.
//! - **A literal** is the shortest decimal that round-trips (`{:e}`) of the series term rounded
//!   once, in integers, at binary64 (`Series::<f64>`): Rust's own parser reads it back to the same
//!   bits, so it is the correctly rounded value of the exact rational. Binary32 waits for its sweep.
//! - **The layout** is what `rustfmt` leaves alone, `just lint` running `cargo fmt --check`, at
//!   its defaults: an array of at most `array_width` (60) characters, brackets excluded, is one line;
//!   a longer one whose elements are all at most `short_array_element_width_threshold` (10)
//!   characters is packed, as many to a line as fit in `max_width` (100); any other is one element
//!   per line. `a_series_is_laid_out_as_rustfmt_leaves_it` runs `rustfmt` over the boundaries.

use sha2::{Digest, Sha256};

use super::{CSV, HEADER};
use crate::seeded::{Coeff, Series, SERIES_FILE};

/// The generated file, relative to the repository root.
pub(super) const PATH: &str = "xtask/src/seeded/generated.rs";

/// rustfmt's default `array_width`: the widest array literal, brackets excluded, on one line.
const ARRAY_WIDTH: usize = 60;
/// rustfmt's default `short_array_element_width_threshold`: an element at most this wide is short.
const SHORT_ELEMENT: usize = 10;
/// rustfmt's default `max_width`: the widest line, trailing comma included.
const MAX_WIDTH: usize = 100;

/// The lower-case hexadecimal SHA-256 of `data`.
fn sha256_hex(data: &[u8]) -> String {
    Sha256::digest(data)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

/// The cell of `row` in the column `name` of [`HEADER`].
fn cell<'a>(row: &'a str, name: &str) -> Result<&'a str, String> {
    let at = HEADER.split(',').position(|c| c == name);
    at.and_then(|i| row.split(',').nth(i))
        .ok_or_else(|| format!("no `{name}` in `{row}`"))
}

/// A finite `f64` written in `{:e}`, as the CSV writes them; the text itself is emitted.
fn number<'a>(row: &'a str, name: &str) -> Result<&'a str, String> {
    let text = cell(row, name)?;
    match text.parse::<f64>() {
        Ok(x) if x.is_finite() => Ok(text),
        _ => Err(format!("{name} `{text}` is not a finite number")),
    }
}

/// `0x3ff0_0000_0000_0000`.
fn grouped(bits: u64) -> String {
    let word = |shift: u32| (bits >> shift) & 0xffff;
    format!(
        "0x{:04x}_{:04x}_{:04x}_{:04x}",
        word(48),
        word(32),
        word(16),
        word(0)
    )
}

/// The series of `terms` literals, at `indent` spaces, as rustfmt lays an array out.
fn series_literal(terms: &[f64], indent: usize) -> String {
    let items: Vec<String> = terms.iter().map(|x| format!("{x:e}")).collect();
    let contents = items.join(", ");
    if contents.len() <= ARRAY_WIDTH {
        return format!("series: [{contents}],");
    }
    let pad = " ".repeat(indent + 4);
    let lines = if items.iter().all(|x| x.len() <= SHORT_ELEMENT) {
        let mut packed: Vec<String> = Vec::new();
        for item in items {
            // The line with this element and the comma that ends it must fit in `max_width`.
            let fits = |line: &String| {
                pad.len() + line.len() + ", ".len() + item.len() + ",".len() <= MAX_WIDTH
            };
            match packed.last_mut() {
                Some(line) if fits(line) => {
                    line.push_str(", ");
                    line.push_str(&item);
                }
                _ => packed.push(item),
            }
        }
        packed
    } else {
        items
    };
    let body: String = lines.iter().map(|line| format!("{pad}{line},\n")).collect();
    format!("series: [\n{body}{}],", " ".repeat(indent))
}

/// The constant of `c` from its CSV `row`.
fn constant(c: Coeff, row: &str, series: &Series<f64>) -> Result<String, String> {
    let name = c.name();
    if cell(row, "coeff")? != name || cell(row, "precision")? != "f64" {
        return Err(format!("`{row}` is not the f64 row of `{name}`"));
    }
    let terms: usize = cell(row, "terms")?
        .parse()
        .map_err(|e| format!("{name}: terms: {e}"))?;
    if !(1..=series.terms()).contains(&terms) {
        return Err(format!(
            "{name}: {terms} terms of a {}-term series",
            series.terms()
        ));
    }
    let bits = cell(row, "switch_bits")?
        .strip_prefix("0x")
        .and_then(|h| u64::from_str_radix(h, 16).ok())
        .ok_or_else(|| format!("{name}: switch_bits is not 0x and 16 hex digits"))?;
    let switch = f64::from_bits(bits);
    // The comment is the bits' own decimal: a CSV whose columns disagree is refused, not emitted.
    let decimal = format!("{switch:e}");
    if !switch.is_finite() || switch < 0.0 || decimal != cell(row, "switch_z")? {
        return Err(format!(
            "{name}: switch_bits and switch_z disagree or are not a z >= 0"
        ));
    }
    let (value, deriv) = (number(row, "value_max_u")?, number(row, "deriv_max_u")?);
    let literals = series_literal(&series.of(c)[..terms], 4);
    Ok(format!(
        "// Objective (max u): value {value}, derivative {deriv}.\n\
         pub(crate) const {upper}_F64: Switch<f64, {terms}> = Switch {{\n    \
         below: f64::from_bits({}), // θ² < {decimal}\n    \
         {literals}\n}};\n",
        grouped(bits),
        upper = name.to_uppercase(),
    ))
}

/// The generated file for the sweep `csv`, with the series of the corpus: `series_file` is the
/// bytes of `coeff_series.jsonl`, of which `series` is the reading.
pub(super) fn render(
    csv: &str,
    series_file: &[u8],
    series: &Series<f64>,
) -> Result<String, String> {
    let mut lines = csv.lines();
    if lines.next() != Some(HEADER) {
        return Err(format!("{CSV} does not start with its documented header"));
    }
    let rows: Vec<&str> = lines.collect();
    if rows.len() != Coeff::ALL.len() {
        return Err(format!(
            "{CSV} has {} rows; expected one per coefficient, {}",
            rows.len(),
            Coeff::ALL.len()
        ));
    }
    let mut out = format!(
        "// @generated by `cargo xtask thresholds` from {CSV} — do not edit.\n\
         // Source sweep rev: sha256 {}.\n\
         // Source series: conformance/corpus/{SERIES_FILE}, sha256 {}.\n\nuse super::switch::Switch;\n",
        sha256_hex(csv.as_bytes()),
        sha256_hex(series_file)
    );
    for (c, row) in Coeff::ALL.into_iter().zip(rows) {
        out.push('\n');
        out.push_str(&constant(c, row, series)?);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::corpus_dir;
    use crate::conformance::root;

    fn committed_csv() -> Result<String, String> {
        let path = root()?.join(CSV);
        std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))
    }

    fn series_file() -> Result<Vec<u8>, String> {
        let path = corpus_dir()?.join(SERIES_FILE);
        std::fs::read(&path).map_err(|e| format!("{}: {e}", path.display()))
    }

    fn series() -> Result<Series<f64>, String> {
        Series::load(&corpus_dir()?)
    }

    /// The switch `synthetic` gives `b`: four different 16-bit words, none zero.
    const B_SWITCH_BITS: u64 = 0x3fc0_1234_5678_9abc;

    /// A sweep CSV of the documented columns with `terms[i]` terms and the switch `2^-(i+1)` for
    /// coefficient `i`, except `b`'s ([`B_SWITCH_BITS`]); the objectives `1.5e0` and `2.5e1`, and
    /// `0` in every other cell.
    fn synthetic(terms: [usize; 6]) -> String {
        let mut out = format!("{HEADER}\n");
        for (i, c) in Coeff::ALL.iter().enumerate() {
            let z = match c {
                Coeff::B => f64::from_bits(B_SWITCH_BITS),
                _ => 1.0 / (1u64 << (i + 1)) as f64,
            };
            let cells: Vec<String> = HEADER
                .split(',')
                .map(|name| match name {
                    "coeff" => c.name().to_string(),
                    "precision" => "f64".to_string(),
                    "terms" => terms[i].to_string(),
                    "switch_bits" => format!("0x{:016x}", z.to_bits()),
                    "switch_z" => format!("{z:e}"),
                    "value_max_u" => "1.5e0".to_string(),
                    "deriv_max_u" => "2.5e1".to_string(),
                    _ => "0".to_string(),
                })
                .collect();
            out.push_str(&cells.join(","));
            out.push('\n');
        }
        out
    }

    #[test]
    fn a_constant_has_the_documented_shape() -> Result<(), String> {
        let text = render(&synthetic([8, 8, 3, 1, 8, 7]), &series_file()?, &series()?)?;
        let mut lines = text.lines();
        let marker = "// @generated by `cargo xtask thresholds` from conformance/sweeps/thresholds.csv — do not edit.";
        assert_eq!(lines.next(), Some(marker));
        for lead in [
            "// Source sweep rev: sha256 ",
            "// Source series: conformance/corpus/coeff_series.jsonl, sha256 ",
        ] {
            let line = lines.next().unwrap_or_default();
            assert!(
                line.starts_with(lead) && line.len() == lead.len() + 64 + 1,
                "{line}"
            );
        }
        assert_eq!(lines.next(), Some(""));
        assert_eq!(lines.next(), Some("use super::switch::Switch;"));
        let lit = |x: f64| format!("{x:e}");
        // `b`: three terms, one literal per line, its switch's four words in order.
        let want = format!(
            "// Objective (max u): value 1.5e0, derivative 2.5e1.\n\
             pub(crate) const B_F64: Switch<f64, 3> = Switch {{\n    \
             below: f64::from_bits(0x3fc0_1234_5678_9abc), // θ² < {}\n    \
             series: [\n        {},\n        {},\n        {},\n    ],\n}};\n",
            lit(f64::from_bits(B_SWITCH_BITS)),
            lit(1.0 / 6.0),
            lit(-1.0 / 120.0),
            lit(1.0 / 5040.0)
        );
        assert!(text.contains(&format!("\n{want}\n")), "{text}");
        // `c`: one term, and a series that fits on one line stays on one.
        let one = "pub(crate) const C_F64: Switch<f64, 1> = Switch {\n    \
            below: f64::from_bits(0x3fb0_0000_0000_0000), // θ² < 6.25e-2\n    \
            series: [8.333333333333333e-2],\n};\n";
        assert!(text.contains(one), "{text}");
        assert_eq!(text.matches("pub(crate) const").count(), 6);
        assert!(text.contains("pub(crate) const E_F64: Switch<f64, 7>") && text.ends_with("};\n"));
        Ok(())
    }

    #[test]
    fn a_switch_is_four_groups_of_sixteen_bits_from_the_high_word_down() {
        assert_eq!(grouped(0x0123_4567_89ab_cdef), "0x0123_4567_89ab_cdef");
        assert_eq!(grouped(0), "0x0000_0000_0000_0000");
        assert_eq!(grouped(1), "0x0000_0000_0000_0001");
        assert_eq!(grouped(u64::MAX), "0xffff_ffff_ffff_ffff");
    }

    #[test]
    fn the_sources_are_named_by_their_digests_and_a_changed_source_changes_the_file(
    ) -> Result<(), String> {
        let (csv, jsonl, series) = (synthetic([8; 6]), series_file()?, series()?);
        let text = render(&csv, &jsonl, &series)?;
        assert_eq!(
            sha256_hex(b"abc"),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        // The series file's digest is the one Python's `hashlib` recorded in the manifest.
        let path = corpus_dir()?.join("MANIFEST.json");
        let manifest = std::fs::read(&path).map_err(|e| format!("{}: {e}", path.display()))?;
        let manifest: serde_json::Value =
            serde_json::from_slice(&manifest).map_err(|e| e.to_string())?;
        let recorded = manifest["files"][SERIES_FILE]["sha256"]
            .as_str()
            .unwrap_or("");
        assert_eq!(sha256_hex(&jsonl), recorded);
        let rev = format!(
            "// Source sweep rev: sha256 {}.\n",
            sha256_hex(csv.as_bytes())
        );
        let source = format!("coeff_series.jsonl, sha256 {recorded}.\n");
        assert!(text.contains(&rev) && text.contains(&source), "{text}");
        // One digit of one objective, or one byte of the series file: the constants are the same,
        // the file is not, and it is the line of that source that moves.
        let other = render(&csv.replacen("1.5e0", "1.6e0", 1), &jsonl, &series)?;
        assert_ne!(other, text);
        assert!(other.contains(&source) && !other.contains(&rev));
        let mut edited = jsonl.clone();
        edited.push(b'\n');
        let other = render(&csv, &edited, &series)?;
        assert!(other.contains(&rev) && !other.contains(&source));
        assert_eq!(render(&csv, &jsonl, &series)?, text);
        Ok(())
    }

    #[test]
    fn every_literal_reads_back_as_the_series_term_rounded_once() -> Result<(), String> {
        let (csv, series) = (committed_csv()?, series()?);
        let text = render(&csv, &series_file()?, &series)?;
        let mut checked = 0;
        let mut expected = 0;
        for (c, row) in Coeff::ALL.into_iter().zip(csv.lines().skip(1)) {
            let terms: usize = cell(row, "terms")?.parse().map_err(|e| format!("{e}"))?;
            let name = format!("pub(crate) const {}_F64", c.name().to_uppercase());
            let body = text.split(&name).nth(1).ok_or("no constant")?;
            let body = body.split("};").next().unwrap_or_default();
            let is_literal = |l: &&str| l.starts_with(|c: char| c == '-' || c.is_ascii_digit());
            let literals: Vec<&str> = body
                .lines()
                .map(str::trim)
                .filter(is_literal)
                .map(|l| l.trim_end_matches(','))
                .collect();
            assert_eq!(literals.len(), terms, "{c:?}");
            for (literal, &want) in literals.iter().zip(series.of(c)) {
                // Rust's own parser, not the harness's integer rounding: the same bits.
                let got: f64 = literal.parse().map_err(|e| format!("{literal}: {e}"))?;
                assert_eq!(got.to_bits(), want.to_bits(), "{c:?}: {literal}");
                checked += 1;
            }
            expected += terms;
        }
        assert_eq!(checked, expected);
        Ok(())
    }

    /// A finite number of which `{:e}` is `width` characters (5 to 19): at most 15 significant
    /// digits, none zero, so the shortest decimal that round-trips is the text itself.
    fn literal(width: usize) -> Result<f64, String> {
        let text = match width {
            5 => "1.5e0".to_string(),
            6..=19 => format!("1.{}e-1", &"12345678901234"[..width - 5]),
            _ => return Err(format!("no literal of width {width}")),
        };
        let x: f64 = text.parse().map_err(|e| format!("{text}: {e}"))?;
        assert_eq!(format!("{x:e}"), text);
        Ok(x)
    }

    /// The constant around a series, as `constant` writes it.
    fn around(terms: &[f64]) -> String {
        format!(
            "pub(crate) const X: Switch<f64, {}> = Switch {{\n    \
             below: f64::from_bits(0x3fc0_0000_0000_0000), // θ² < 1.25e-1\n    {}\n}};\n",
            terms.len(),
            series_literal(terms, 4)
        )
    }

    /// The diff `rustfmt --check` reports for `text`, at the workspace edition and the defaults of
    /// this repository (no `rustfmt.toml`); `None` when it leaves `text` as it is.
    fn rustfmt_diff(text: &str) -> Result<Option<String>, String> {
        use std::io::Write;
        use std::process::{Command, Stdio};
        let mut child = Command::new("rustfmt")
            .args(["--edition", "2021", "--check"])
            .current_dir(root()?)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| format!("rustfmt: {e} (`rustup component add rustfmt`)"))?;
        let stdin = child.stdin.as_mut().ok_or("no stdin")?;
        stdin
            .write_all(text.as_bytes())
            .map_err(|e| e.to_string())?;
        let out = child.wait_with_output().map_err(|e| e.to_string())?;
        // On stdin `--check` prints the diff and exits 0; on a syntax error it prints nothing.
        let diff = String::from_utf8_lossy(&out.stdout);
        match (diff.is_empty(), out.status.success()) {
            (true, true) => Ok(None),
            (true, false) => Err(format!("rustfmt: {}", String::from_utf8_lossy(&out.stderr))),
            (false, _) => Ok(Some(diff.lines().take(40).collect::<Vec<_>>().join("\n"))),
        }
    }

    #[test]
    fn a_series_is_laid_out_as_rustfmt_leaves_it() -> Result<(), String> {
        // Every mix of up to four widths around the 60-character limit and the 10-character
        // threshold, every uniform series up to the corpus's 16 terms, and the packed line's
        // limit, in one file.
        let widths = [5, 10, 11, 18, 19];
        let mut mixes: Vec<Vec<usize>> = vec![Vec::new()];
        let mut all: Vec<Vec<usize>> = Vec::new();
        for _ in 0..4 {
            mixes = mixes
                .iter()
                .flat_map(|m| widths.iter().map(|&w| [m.as_slice(), &[w]].concat()))
                .collect();
            all.extend(mixes.iter().cloned());
        }
        for w in [5, 8, 10, 11, 14, 19] {
            all.extend((1..=16).map(|n| vec![w; n]));
        }
        // Twelve `1.5e0` and one more element of 6 to 10 characters: a line of 99 to 103.
        all.extend((6..=10).map(|w| [vec![5; 12], vec![w]].concat()));
        let contents = |m: &Vec<usize>| m.iter().sum::<usize>() + 2 * (m.len() - 1);
        for edge in [59, 60, 61, 62] {
            assert!(
                all.iter().any(|m| contents(m) == edge),
                "no series of {edge}"
            );
        }
        let short = |m: &&Vec<usize>| m.iter().all(|&w| w <= SHORT_ELEMENT);
        assert!(all.iter().filter(short).any(|m| contents(m) > ARRAY_WIDTH));
        let mut text = String::new();
        for mix in &all {
            let terms: Result<Vec<f64>, String> = mix.iter().map(|&w| literal(w)).collect();
            text.push_str(&around(&terms?));
        }
        assert_eq!(rustfmt_diff(&text)?, None);
        // The check can fail: a series that is one per line where rustfmt packs it.
        let unpacked = around(&[1.5; 9]).replace("1.5e0, ", "1.5e0,\n        ");
        assert!(rustfmt_diff(&unpacked)?.is_some(), "{unpacked}");
        Ok(())
    }

    #[test]
    fn a_short_series_is_one_line_and_a_long_one_is_packed_or_one_per_line() -> Result<(), String> {
        assert_eq!(series_literal(&[0.5], 4), "series: [5e-1],");
        assert_eq!(series_literal(&[0.5, -0.25], 4), "series: [5e-1, -2.5e-1],");
        let long = series_literal(&[1.0 / 3.0; 4], 4);
        assert_eq!(long.lines().count(), 6, "{long}");
        assert!(long.starts_with("series: [\n        3.333333333333333e-1,\n"));
        assert!(long.ends_with("\n    ],"));
        // `array_width` counts the elements and their separators, not the brackets: eight `1.5e0`
        // are 54 characters.
        assert_eq!(series_literal(&[1.5; 8], 4).lines().count(), 1);
        // Elements that are all short are packed once past it, thirteen to a line at `max_width`.
        let packed = series_literal(&[1.5; 16], 4);
        let row = |n: usize| format!("        {},\n", ["1.5e0"; 16][..n].join(", "));
        assert_eq!(packed, format!("series: [\n{}{}    ],", row(13), row(3)));
        // Long elements, 60 characters and one past them; and a single long one among short
        // ones makes the array one per line, not packed.
        let (a, b, c) = (literal(19)?, literal(18)?, literal(5)?);
        assert_eq!(series_literal(&[a, a, b], 4).lines().count(), 1);
        assert_eq!(series_literal(&[a, a, a], 4).lines().count(), 5);
        assert_eq!(
            series_literal(&[c, c, c, c, c, c, c, c, a], 4)
                .lines()
                .count(),
            11
        );
        Ok(())
    }

    #[test]
    fn a_csv_that_is_not_the_documented_one_is_refused() -> Result<(), String> {
        let (jsonl, series) = (series_file()?, series()?);
        let base = synthetic([8; 6]);
        let refuses = |csv: &str, why: &str| {
            let e = render(csv, &jsonl, &series).err().unwrap_or_default();
            assert!(e.contains(why), "{why}: {e}");
        };
        let edit = |from: &str, to: &str| {
            assert!(base.contains(from), "{from}");
            base.replacen(from, to, 1)
        };
        refuses("", "documented header");
        refuses(&base.replacen("coeff,", "coeff ,", 1), "documented header");
        let first = base.lines().take(2).collect::<Vec<_>>().join("\n") + "\n";
        refuses(&first, "expected one per coefficient, 6");
        refuses(&edit("c,f64,8", "c,f32,8"), "f64 row of `c`");
        refuses(&edit("k,f64,8", "c,f64,8"), "f64 row of `k`");
        refuses(&edit("c,f64,8", "c,f64,17"), "17 terms");
        refuses(&edit("c,f64,8", "c,f64,0"), "0 terms");
        refuses(&edit("c,f64,8", "c,f64,x"), "terms");
        // 2^-2 is 0x3fd0..; the decimal column must be the bits' own, and a z below 0 is refused.
        refuses(
            &edit("0x3fd0000000000000", "0x3fd0000000000001"),
            "disagree",
        );
        refuses(
            &edit("0x3fd0000000000000", "3fd0000000000000"),
            "16 hex digits",
        );
        refuses(
            &edit("0x3fd0000000000000", "0xbfd0000000000000"),
            "disagree",
        );
        refuses(&edit("1.5e0", "inf"), "not a finite number");
        refuses(&edit("2.5e1", "x"), "not a finite number");
        Ok(())
    }
}
