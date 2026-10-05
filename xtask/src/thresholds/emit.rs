//! The generated files of `docs/PHASE1.md` §6 (`0004` item 3, `0016` item 3): the seeded kernels'
//! `xtask/src/seeded/generated.rs` and `helicoid::coeffs`'
//! `crates/helicoid/src/coeffs/generated.rs`, each a function of its sweep CSV and
//! `coeff_series.jsonl` and of nothing else, so `thresholds --check` compares them byte for byte.
//! One `Switch` per coefficient and precision in the CSV's order, `K_F64` to `R_F64`, then `K_F32` to
//! `R_F32`. The seeded file's hold their series (`seeded::switch`); `helicoid`'s take theirs from
//! the swept series that follow them in the same file, `SWEPT_K_F64` to `SWEPT_R_F32`, the first
//! `TERMS` terms of each row (`coeffs::Switch::first`), which is what the sweep of that target
//! measures, so the file is a fixed point (`render_helicoid`).
//!
//! Each file is compiled into this tool (`seeded:correct` runs the first, the sweep measures the
//! second's swept series), so a hand edit that no longer compiles stops `cargo xtask` from building,
//! and neither `thresholds` nor `thresholds --check` can name it: `git restore` the file first. One
//! that compiles is named, file and line.
//!
//! Readings where §6 is silent, each the smallest (the first two are 0014 (draft) questions 12 and 13):
//!
//! - **The source revision** is two SHA-256 digests, of the CSV's bytes and of the series file's,
//!   which supplies the literals. A git revision cannot name the commit that contains the file and
//!   differs on every commit, which no byte-for-byte check survives.
//! - **The objective line** is one per constant, above it, from the CSV's own columns: §6's
//!   example shows one line for one constant, and a file holds sixteen.
//! - **`terms`** is the array length, so the kernel's `terms` is `series.len()` and the constant
//!   holds the terms the sweep measured and no others.
//! - **A literal** is the shortest decimal that round-trips (`{:e}`) of the series term rounded
//!   once, in integers, at the constant's own precision (`Series::<f64>`, `Series::<f32>`): Rust's
//!   own parser reads it back to the same bits, so it is the correctly rounded value of the exact
//!   rational, and a binary32 literal is never a binary64 one rounded again (`0016` item 3).
//! - **`r`'s constant** is the series in `s = n²/w²` (the comment names it); its prefactor `2/w`
//!   is the kernel's.
//! - **The layout** is what `rustfmt` leaves alone, `just lint` running `cargo fmt --check`, at
//!   its defaults: an array of at most `array_width` (60) characters, brackets excluded, is one line;
//!   a longer one whose elements are all at most `short_array_element_width_threshold` (10)
//!   characters is packed, as many to a line as fit in `max_width` (100); any other is one element
//!   per line. `a_series_is_laid_out_as_rustfmt_leaves_it` runs `rustfmt` over the boundaries.

use std::fmt::LowerExp;

use helicoid_linalg::Precision;
use sha2::{Digest, Sha256};

use super::search::TERMS;
use super::{precision_name, shown, CSV_HELICOID, CSV_SEEDED, HEADER};
use crate::seeded::{Series, Swept, SERIES_FILE};

/// The seeded kernels' generated file, relative to the repository root.
pub(super) const PATH_SEEDED: &str = "xtask/src/seeded/generated.rs";
/// `helicoid::coeffs`' generated file, relative to the repository root.
pub(super) const PATH_HELICOID: &str = "crates/helicoid/src/coeffs/generated.rs";

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

/// `bits` as `hex` hex digits in groups of four: `0x3ff0_0000_0000_0000`, `0x3f80_0000`.
fn grouped(bits: u64, hex: usize) -> String {
    let digits = format!("{bits:0hex$x}");
    let words: Vec<&str> = (0..hex).step_by(4).map(|i| &digits[i..i + 4]).collect();
    format!("0x{}", words.join("_"))
}

/// The array `head` `terms` `tail` (`series: [` ... `],`, `const X: [f64; 8] = [` ... `];`) at
/// `indent` spaces, as rustfmt lays an array out.
fn array_literal<T: LowerExp>(head: &str, terms: &[T], indent: usize, tail: &str) -> String {
    let items: Vec<String> = terms.iter().map(|x| format!("{x:e}")).collect();
    let contents = items.join(", ");
    let one_line = indent + head.chars().count() + contents.len() + tail.len();
    if contents.len() <= ARRAY_WIDTH {
        if one_line <= MAX_WIDTH {
            return format!("{head}{contents}{tail}");
        }
        // A constant whose head is too long for one line: rustfmt breaks after the `=`.
        if let Some(lhs) = head.strip_suffix('[') {
            let pad = " ".repeat(indent + 4);
            return format!("{}\n{pad}[{contents}{}", lhs.trim_end(), tail);
        }
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
    format!("{head}\n{body}{}{tail}", " ".repeat(indent))
}

/// The series of `terms` literals in a `Switch`, at `indent` spaces.
fn series_literal<T: LowerExp>(terms: &[T], indent: usize) -> String {
    array_literal("series: [", terms, indent, "],")
}

/// What a CSV row chose for a coefficient at a precision, checked against the row's own columns.
struct Chosen<'a> {
    terms: usize,
    /// The bit pattern, with the `hex` digits of the precision.
    bits: u64,
    hex: usize,
    /// The switch's decimal at its precision, as the comment prints it.
    decimal: String,
    value: &'a str,
    deriv: &'a str,
    /// The search space the choice was made in, as one comment line (`docs/decisions/0039`
    /// item 9). An objective reported alone reads as an optimum even when it is the edge of the
    /// box, which is how eight switches sat on a wall for the file's whole existence.
    feasibility: String,
}

/// The choice in `row`, the `precision` row of `id`, of at most `max_terms` terms.
fn chosen<'a>(
    id: Swept,
    precision: Precision,
    row: &'a str,
    max_terms: usize,
) -> Result<Chosen<'a>, String> {
    let (name, p) = (id.name(), precision_name(precision));
    if cell(row, "coeff")? != name || cell(row, "precision")? != p {
        return Err(format!("`{row}` is not the {p} row of `{name}`"));
    }
    let terms: usize = cell(row, "terms")?
        .parse()
        .map_err(|e| format!("{name}: terms: {e}"))?;
    if !(1..=max_terms).contains(&terms) {
        return Err(format!(
            "{name}: {terms} terms of a {max_terms}-term series"
        ));
    }
    let hex = match precision {
        Precision::F64 => 16,
        Precision::F32 => 8,
    };
    let bits = cell(row, "switch_bits")?
        .strip_prefix("0x")
        .filter(|h| h.len() == hex)
        .and_then(|h| u64::from_str_radix(h, 16).ok())
        .ok_or_else(|| format!("{name}: switch_bits is not 0x and {hex} hex digits"))?;
    let switch = match precision {
        Precision::F64 => f64::from_bits(bits),
        Precision::F32 => f64::from(f32::from_bits(bits as u32)),
    };
    // The comment is the bits' own decimal: a CSV whose columns disagree is refused, not emitted.
    let decimal = shown(precision, switch);
    if !switch.is_finite() || switch < 0.0 || decimal != cell(row, "switch_z")? {
        return Err(format!(
            "{name}: switch_bits and switch_z disagree or are not a z >= 0"
        ));
    }
    let binding = cell(row, "binding")?;
    if !["none", "domain", "terms", "domain+terms"].contains(&binding) {
        return Err(format!(
            "{name}: binding `{binding}` is not one of the four"
        ));
    }
    Ok(Chosen {
        terms,
        bits,
        hex,
        decimal,
        value: number(row, "value_max_u")?,
        deriv: number(row, "deriv_max_u")?,
        feasibility: format!(
            "grid [{}, {}] of {} per decade searched below π², terms <= {}; binding: {binding}",
            cell(row, "grid_lo")?,
            cell(row, "grid_hi")?,
            64,
            cell(row, "term_cap")?,
        ),
    })
}

/// The name of the constant of `id` at `precision`, `K_F64`.
fn const_name(id: Swept, precision: Precision) -> String {
    let p = precision_name(precision);
    format!("{}_{}", id.name().to_uppercase(), p.to_uppercase())
}

/// The variable a switch is in.
fn variable(id: Swept) -> &'static str {
    if id == Swept::R {
        "n²/w²"
    } else {
        "θ²"
    }
}

/// The seeded constant of `id` at `precision` from its CSV `row`.
fn constant(
    id: Swept,
    precision: Precision,
    row: &str,
    (wide, narrow): (&Series<f64>, &Series<f32>),
) -> Result<String, String> {
    let c = chosen(id, precision, row, wide.terms())?;
    let literals = match precision {
        Precision::F64 => series_literal(&wide.swept(id)[..c.terms], 4),
        Precision::F32 => series_literal(&narrow.swept(id)[..c.terms], 4),
    };
    let p = precision_name(precision);
    Ok(format!(
        "// Objective (max u): value {}, derivative {}.\n\
         // Feasibility: {}.\n\
         pub(crate) const {}: Switch<{p}, {}> = Switch {{\n    \
         below: {p}::from_bits({}), // {} < {}\n    \
         {literals}\n}};\n",
        c.value,
        c.deriv,
        c.feasibility,
        const_name(id, precision),
        c.terms,
        grouped(c.bits, c.hex),
        variable(id),
        c.decimal,
    ))
}

/// The first three lines of a generated file: the marker and the two digests it is a function of.
fn header(csv_path: &str, csv: &str, series_file: &[u8]) -> String {
    format!(
        "// @generated by `cargo xtask thresholds` from {csv_path} — do not edit.\n\
         // Source sweep rev: sha256 {}.\n\
         // Source series: conformance/corpus/{SERIES_FILE}, sha256 {}.\n",
        sha256_hex(csv.as_bytes()),
        sha256_hex(series_file)
    )
}

/// The data rows of `csv`, the sweep at `csv_path`: its header is the documented one and it has a
/// row per coefficient and precision.
fn rows<'a>(csv_path: &str, csv: &'a str) -> Result<Vec<&'a str>, String> {
    let mut lines = csv.lines();
    if lines.next() != Some(HEADER) {
        return Err(format!(
            "{csv_path} does not start with its documented header"
        ));
    }
    let rows: Vec<&str> = lines.collect();
    let of = Swept::ALL.len();
    if rows.len() != 2 * of {
        return Err(format!(
            "{csv_path} has {} rows; expected one per coefficient and precision, {}",
            rows.len(),
            2 * of
        ));
    }
    Ok(rows)
}

/// The precision and coefficient of row `i` of the CSV.
fn row_of(i: usize) -> (Precision, Swept) {
    let of = Swept::ALL.len();
    ([Precision::F64, Precision::F32][i / of], Swept::ALL[i % of])
}

/// The seeded generated file for the sweep `csv` (at [`CSV_SEEDED`]), with the series of the
/// corpus at both precisions: `series_file` is the bytes of `coeff_series.jsonl`, of which the
/// series are the readings.
pub(super) fn render(
    csv: &str,
    series_file: &[u8],
    wide: &Series<f64>,
    narrow: &Series<f32>,
) -> Result<String, String> {
    let rows = rows(CSV_SEEDED, csv)?;
    let mut out = header(CSV_SEEDED, csv, series_file);
    out.push_str("\nuse super::switch::Switch;\n");
    for (i, row) in rows.into_iter().enumerate() {
        let (precision, id) = row_of(i);
        out.push('\n');
        out.push_str(&constant(id, precision, row, (wide, narrow))?);
    }
    Ok(out)
}

/// The generated file of `helicoid::coeffs` for the sweep `csv` (at [`CSV_HELICOID`]): a `Switch`
/// per coefficient and precision, then the swept series they take their terms from. `placeholder`
/// marks a file whose switches are not a sweep's (`bootstrap`), written when the compiled-in swept
/// series are not the corpus's.
pub(super) fn render_helicoid(
    csv: &str,
    series_file: &[u8],
    (wide, narrow): (&Series<f64>, &Series<f32>),
    placeholder: bool,
) -> Result<String, String> {
    let rows = rows(CSV_HELICOID, csv)?;
    let mut out = header(CSV_HELICOID, csv, series_file);
    if placeholder {
        out.push_str(
            "// Placeholder switches: the series changed; run `cargo xtask thresholds` again.\n",
        );
    }
    out.push_str("\nuse super::Switch;\n");
    for (i, row) in rows.into_iter().enumerate() {
        let (precision, id) = row_of(i);
        let c = chosen(id, precision, row, TERMS)?;
        let p = precision_name(precision);
        out.push_str(&format!(
            "\n// Objective (max u): value {}, derivative {}.\n\
             // Feasibility: {}.\n\
             pub(crate) const {name}: Switch<{p}, {}> = Switch::first(\n    \
             {p}::from_bits({}), // {} < {}\n    \
             &SWEPT_{name},\n);\n",
            c.value,
            c.deriv,
            c.feasibility,
            c.terms,
            grouped(c.bits, c.hex),
            variable(id),
            c.decimal,
            name = const_name(id, precision),
        ));
    }
    out.push_str(&format!(
        "\n// The series `helicoid::__sweep` measures: the first {TERMS} terms of each row, the exact\n\
         // rational rounded once at its precision. A function of the corpus alone.\n"
    ));
    for i in 0..2 * Swept::ALL.len() {
        let (precision, id) = row_of(i);
        let p = precision_name(precision);
        let name = const_name(id, precision);
        let head = format!("pub(crate) const SWEPT_{name}: [{p}; {TERMS}] = [");
        let literal = match precision {
            Precision::F64 => array_literal(&head, &wide.swept(id)[..TERMS], 0, "];"),
            Precision::F32 => array_literal(&head, &narrow.swept(id)[..TERMS], 0, "];"),
        };
        out.push_str(&format!("{}{literal}\n", if i == 0 { "" } else { "\n" }));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::corpus_dir;
    use crate::conformance::number::Dyadic;
    use crate::conformance::root;
    use crate::conformance::testkit::Scratch;
    use num_bigint::BigUint;

    fn committed_csv() -> Result<String, String> {
        let path = root()?.join(CSV_SEEDED);
        std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))
    }

    fn series_file() -> Result<Vec<u8>, String> {
        let path = corpus_dir()?.join(SERIES_FILE);
        std::fs::read(&path).map_err(|e| format!("{}: {e}", path.display()))
    }

    fn series() -> Result<(Series<f64>, Series<f32>), String> {
        let dir = corpus_dir()?;
        Ok((Series::load(&dir)?, Series::load(&dir)?))
    }

    /// The switch `synthetic` gives `b` at `f64`: four different 16-bit words, none zero.
    const B_SWITCH_BITS: u64 = 0x3fc0_1234_5678_9abc;
    /// And at `f32`: two.
    const B_SWITCH_BITS32: u32 = 0x3c01_2345;

    /// A sweep CSV of the documented columns with `terms[i]` terms and the switch `2^-(i+1)` for
    /// each id `i` at both precisions, except `b`'s (`B_SWITCH_BITS`, `B_SWITCH_BITS32`); the
    /// objectives `1.5e0` and `2.5e1`, and `0` in every other cell.
    fn synthetic(terms: [usize; 8]) -> String {
        let mut out = format!("{HEADER}\n");
        for precision in [Precision::F64, Precision::F32] {
            for (i, id) in Swept::ALL.iter().enumerate() {
                let z = 1.0 / (1u64 << (i + 1)) as f64;
                let (bits, z) = match (precision, id) {
                    (Precision::F64, Swept::Coeff(crate::seeded::Coeff::B)) => {
                        let z = f64::from_bits(B_SWITCH_BITS);
                        (format!("0x{:016x}", B_SWITCH_BITS), z)
                    }
                    (Precision::F32, Swept::Coeff(crate::seeded::Coeff::B)) => {
                        let z = f64::from(f32::from_bits(B_SWITCH_BITS32));
                        (format!("0x{:08x}", B_SWITCH_BITS32), z)
                    }
                    (Precision::F64, _) => (format!("0x{:016x}", z.to_bits()), z),
                    (Precision::F32, _) => (format!("0x{:08x}", (z as f32).to_bits()), z),
                };
                let cells: Vec<String> = HEADER
                    .split(',')
                    .map(|name| match name {
                        "coeff" => id.name().to_string(),
                        "precision" => precision_name(precision).to_string(),
                        "terms" => terms[i].to_string(),
                        "switch_bits" => bits.clone(),
                        "switch_z" => shown(precision, z),
                        "value_max_u" => "1.5e0".to_string(),
                        "deriv_max_u" => "2.5e1".to_string(),
                        "grid_lo" => "1e-16".to_string(),
                        "grid_hi" => "1e1".to_string(),
                        "term_cap" => TERMS.to_string(),
                        "binding" => "none".to_string(),
                        _ => "0".to_string(),
                    })
                    .collect();
                out.push_str(&cells.join(","));
                out.push('\n');
            }
        }
        out
    }

    /// `x` at `f32`, as the series holds it.
    fn lit32(x: f32) -> String {
        format!("{x:e}")
    }

    #[test]
    fn a_constant_has_the_documented_shape() -> Result<(), String> {
        let (wide, narrow) = series()?;
        let text = render(
            &synthetic([8, 8, 3, 1, 8, 7, 8, 5]),
            &series_file()?,
            &wide,
            &narrow,
        )?;
        let mut lines = text.lines();
        let marker = "// @generated by `cargo xtask thresholds` from conformance/sweeps/thresholds-seeded.csv — do not edit.";
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
             // Feasibility: grid [1e-16, 1e1] of 64 per decade searched below π², terms <= 16; \
             binding: none.\n\
             pub(crate) const B_F64: Switch<f64, 3> = Switch {{\n    \
             below: f64::from_bits(0x3fc0_1234_5678_9abc), // θ² < {}\n    \
             series: [\n        {},\n        {},\n        {},\n    ],\n}};\n",
            lit(f64::from_bits(B_SWITCH_BITS)),
            lit(1.0 / 6.0),
            lit(-1.0 / 120.0),
            lit(1.0 / 5040.0)
        );
        assert!(text.contains(&format!("\n{want}\n")), "{text}");
        // The same at `f32`: two words of bits, literals of binary32's own rounding.
        let want32 = format!(
            "// Objective (max u): value 1.5e0, derivative 2.5e1.\n\
             // Feasibility: grid [1e-16, 1e1] of 64 per decade searched below π², terms <= 16; \
             binding: none.\n\
             pub(crate) const B_F32: Switch<f32, 3> = Switch {{\n    \
             below: f32::from_bits(0x3c01_2345), // θ² < {}\n    \
             series: [{}, {}, {}],\n}};\n",
            lit32(f32::from_bits(B_SWITCH_BITS32)),
            lit32(1.0 / 6.0),
            lit32(-1.0 / 120.0),
            lit32(1.0 / 5040.0)
        );
        assert!(text.contains(&format!("\n{want32}\n")), "{text}");
        // `c`: one term, and a series that fits on one line stays on one.
        let one = "pub(crate) const C_F64: Switch<f64, 1> = Switch {\n    \
            below: f64::from_bits(0x3fb0_0000_0000_0000), // θ² < 6.25e-2\n    \
            series: [8.333333333333333e-2],\n};\n";
        assert!(text.contains(one), "{text}");
        // `r` is a series in `n²/w²`, and its name is one letter.
        assert!(
            text.contains("pub(crate) const R_F32: Switch<f32, 5>"),
            "{text}"
        );
        assert!(text.contains(", // n²/w² < 3.90625e-3\n"), "{text}");
        assert_eq!(text.matches("pub(crate) const").count(), 16);
        assert!(text.contains("pub(crate) const COS_HALF_F64: Switch<f64, 8>"));
        assert!(text.contains("pub(crate) const E_F64: Switch<f64, 7>") && text.ends_with("};\n"));
        Ok(())
    }

    #[test]
    fn a_switch_is_groups_of_sixteen_bits_from_the_high_word_down() {
        assert_eq!(grouped(0x0123_4567_89ab_cdef, 16), "0x0123_4567_89ab_cdef");
        assert_eq!(grouped(0, 16), "0x0000_0000_0000_0000");
        assert_eq!(grouped(1, 16), "0x0000_0000_0000_0001");
        assert_eq!(grouped(u64::MAX, 16), "0xffff_ffff_ffff_ffff");
        assert_eq!(grouped(0x3f80_0000, 8), "0x3f80_0000");
        assert_eq!(grouped(1, 8), "0x0000_0001");
    }

    #[test]
    fn the_sources_are_named_by_their_digests_and_a_changed_source_changes_the_file(
    ) -> Result<(), String> {
        let (csv, jsonl, (wide, narrow)) = (synthetic([8; 8]), series_file()?, series()?);
        let text = render(&csv, &jsonl, &wide, &narrow)?;
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
        let other = render(&csv.replacen("1.5e0", "1.6e0", 1), &jsonl, &wide, &narrow)?;
        assert_ne!(other, text);
        assert!(other.contains(&source) && !other.contains(&rev));
        let mut edited = jsonl.clone();
        edited.push(b'\n');
        let other = render(&csv, &edited, &wide, &narrow)?;
        assert!(other.contains(&rev) && !other.contains(&source));
        assert_eq!(render(&csv, &jsonl, &wide, &narrow)?, text);
        Ok(())
    }

    /// The literals of the constant `name` in `text`, one line or one per line.
    fn literals_of<'a>(text: &'a str, name: &str) -> Result<Vec<&'a str>, String> {
        let body = text.split(&format!("pub(crate) const {name}:")).nth(1);
        let body = body
            .ok_or("no constant")?
            .split("};")
            .next()
            .unwrap_or_default();
        let list = body.split("series: [").nth(1).ok_or("no series")?;
        let list = list.split(']').next().unwrap_or_default();
        Ok(list
            .split(',')
            .map(str::trim)
            .filter(|x| !x.is_empty())
            .collect())
    }

    #[test]
    fn every_literal_reads_back_as_the_series_term_rounded_once_at_its_precision(
    ) -> Result<(), String> {
        let (csv, (wide, narrow)) = (committed_csv()?, series()?);
        let text = render(&csv, &series_file()?, &wide, &narrow)?;
        let (mut checked, mut expected) = (0, 0);
        for (i, row) in csv.lines().skip(1).enumerate() {
            let id = Swept::ALL[i % 8];
            let terms: usize = cell(row, "terms")?.parse().map_err(|e| format!("{e}"))?;
            let p = cell(row, "precision")?;
            let literals = literals_of(
                &text,
                &format!("{}_{}", id.name().to_uppercase(), p.to_uppercase()),
            )?;
            assert_eq!(literals.len(), terms, "{id:?} {p}");
            for (j, literal) in literals.iter().enumerate() {
                // Rust's own parser, not the harness's integer rounding: the same bits, at the
                // precision of the constant.
                let (got, want) = match p {
                    "f64" => (literal.parse::<f64>(), wide.swept(id)[j]),
                    _ => (
                        literal.parse::<f32>().map(f64::from),
                        f64::from(narrow.swept(id)[j]),
                    ),
                };
                let got = got.map_err(|e| format!("{literal}: {e}"))?;
                assert_eq!(got.to_bits(), want.to_bits(), "{id:?} {p}: {literal}");
                checked += 1;
            }
            expected += terms;
        }
        assert_eq!((checked, csv.lines().count()), (expected, 17));
        Ok(())
    }

    /// `x`, a finite binary64, times `2^1074`: an integer.
    fn scaled(x: f64) -> BigUint {
        let d = Dyadic::of(x);
        BigUint::from(d.mant) << (d.exp + 1074) as usize
    }

    /// `literal` is the binary32 nearest to the exact rational `term` (`"num/den"`). Not through
    /// `ratio_to_f32`: `num/den` lies between the midpoints from the literal to its two binary32
    /// neighbours, in integers (`2 num 2^1074` against `den (a + b)`, `a`, `b` the values times
    /// `2^1074`), so no rounding of a binary64 value could pass.
    fn is_the_nearest_binary32(literal: &str, term: &serde_json::Value) -> Result<bool, String> {
        let term = term.as_str().and_then(|t| t.split_once('/'));
        let big = |t: &str| BigUint::parse_bytes(t.trim_start_matches('-').as_bytes(), 10);
        let (num, den) = term
            .and_then(|(n, d)| Some((big(n)?, big(d)?)))
            .ok_or("bad term")?;
        let x = literal
            .parse::<f32>()
            .map_err(|e| format!("{literal}: {e}"))?
            .abs();
        let (lo, hi) = (
            scaled(f64::from(x.next_down())),
            scaled(f64::from(x.next_up())),
        );
        let (here, twice) = (scaled(f64::from(x)), (num << 1075usize));
        Ok(&den * (&lo + &here) <= twice && twice <= &den * (&here + &hi))
    }

    #[test]
    fn every_f32_literal_is_the_nearest_binary32_to_its_exact_rational() -> Result<(), String> {
        let (csv, file) = (committed_csv()?, series_file()?);
        let text = render(&csv, &file, &series()?.0, &series()?.1)?;
        let file = String::from_utf8(file).map_err(|e| e.to_string())?;
        let mut checked = 0;
        for line in file.lines() {
            let row: serde_json::Value = serde_json::from_str(line).map_err(|e| e.to_string())?;
            let name = row["coeff"].as_str().ok_or("no coeff")?;
            let literals = literals_of(&text, &format!("{}_F32", name.to_uppercase()))?;
            for (term, literal) in row["series"]
                .as_array()
                .ok_or("no series")?
                .iter()
                .zip(literals)
            {
                assert!(is_the_nearest_binary32(literal, term)?, "{name}: {literal}");
                checked += 1;
            }
        }
        let mut f32_terms = 0;
        for row in csv.lines().skip(9) {
            f32_terms += cell(row, "terms")?
                .parse::<usize>()
                .map_err(|e| e.to_string())?;
        }
        assert_eq!(checked, f32_terms);
        Ok(())
    }

    /// The literals of the swept series `SWEPT_<name>` in a `helicoid` file.
    fn swept_of<'a>(text: &'a str, name: &str) -> Result<Vec<&'a str>, String> {
        let body = text
            .split(&format!("pub(crate) const SWEPT_{name}: ["))
            .nth(1);
        let list = body
            .and_then(|b| b.split("= [").nth(1))
            .ok_or("no series")?;
        let list = list.split("];").next().unwrap_or_default();
        Ok(list
            .split(',')
            .map(str::trim)
            .filter(|x| !x.is_empty())
            .collect())
    }

    #[test]
    fn a_helicoid_file_has_the_switches_and_the_series_they_take_their_terms_from(
    ) -> Result<(), String> {
        let (wide, narrow) = series()?;
        let jsonl = series_file()?;
        let csv = synthetic([8, 8, 3, 1, 8, 7, 8, 5]);
        let text = render_helicoid(&csv, &jsonl, (&wide, &narrow), false)?;
        let mut lines = text.lines();
        let marker = "// @generated by `cargo xtask thresholds` from conformance/sweeps/thresholds.csv — do not edit.";
        assert_eq!(lines.next(), Some(marker));
        assert!(lines
            .next()
            .is_some_and(|l| l.starts_with("// Source sweep rev: sha256 ")));
        assert!(lines
            .next()
            .is_some_and(|l| l.starts_with("// Source series: ")));
        assert_eq!(lines.next(), Some(""));
        assert_eq!(lines.next(), Some("use super::Switch;"));
        // `b`: three terms of the swept eight, its switch's words in order, at both precisions.
        let want = format!(
            "// Objective (max u): value 1.5e0, derivative 2.5e1.\n\
             // Feasibility: grid [1e-16, 1e1] of 64 per decade searched below π², terms <= 16; \
             binding: none.\n\
             pub(crate) const B_F64: Switch<f64, 3> = Switch::first(\n    \
             f64::from_bits(0x3fc0_1234_5678_9abc), // θ² < {:e}\n    \
             &SWEPT_B_F64,\n);\n",
            f64::from_bits(B_SWITCH_BITS)
        );
        assert!(text.contains(&format!("\n{want}\n")), "{text}");
        assert!(text.contains(&format!(
            "f32::from_bits(0x3c01_2345), // θ² < {}\n    &SWEPT_B_F32,\n);\n",
            lit32(f32::from_bits(B_SWITCH_BITS32))
        )));
        // `r`'s is in `n²/w²`, and every constant is there once: 16 switches, then 16 series.
        assert!(text.contains(", // n²/w² < 3.90625e-3\n"), "{text}");
        assert_eq!(text.matches("pub(crate) const").count(), 32);
        assert_eq!(text.matches("Switch::first(").count(), 16);
        let switches = text
            .find("pub(crate) const SWEPT_K_F64")
            .ok_or("no series")?;
        assert!(text[..switches].contains("_F32: Switch<f32,") && text.ends_with("];\n"));
        // A file that is not a sweep's says so, and a series longer than the swept ones is refused.
        assert!(!text.contains("Placeholder"));
        let holder = render_helicoid(&csv, &jsonl, (&wide, &narrow), true)?;
        assert!(holder
            .lines()
            .nth(3)
            .is_some_and(|l| l.starts_with("// Placeholder switches")));
        let over = render_helicoid(&synthetic([TERMS + 1; 8]), &jsonl, (&wide, &narrow), false);
        let want = format!("{} terms of a {TERMS}-term series", TERMS + 1);
        assert!(over.err().is_some_and(|e| e.contains(&want)));
        Ok(())
    }

    #[test]
    fn every_swept_literal_reads_back_as_the_series_term_rounded_once_at_its_precision(
    ) -> Result<(), String> {
        let (wide, narrow) = series()?;
        let file = String::from_utf8(series_file()?).map_err(|e| e.to_string())?;
        let text = render_helicoid(&committed_csv()?, file.as_bytes(), (&wide, &narrow), false)?;
        let mut checked = 0;
        for line in file.lines() {
            let row: serde_json::Value = serde_json::from_str(line).map_err(|e| e.to_string())?;
            let name = row["coeff"].as_str().ok_or("no coeff")?;
            let id = Swept::ALL
                .into_iter()
                .find(|s| s.name() == name)
                .ok_or("unknown coeff")?;
            let (w, n) = (
                swept_of(&text, &format!("{}_F64", name.to_uppercase()))?,
                swept_of(&text, &format!("{}_F32", name.to_uppercase()))?,
            );
            assert_eq!((w.len(), n.len()), (TERMS, TERMS), "{name}");
            for j in 0..TERMS {
                let got = (w[j].parse::<f64>(), n[j].parse::<f32>());
                let got = (
                    got.0.map_err(|e| e.to_string())?,
                    got.1.map_err(|e| e.to_string())?,
                );
                assert_eq!(got.0.to_bits(), wide.swept(id)[j].to_bits(), "{name}[{j}]");
                assert_eq!(
                    got.1.to_bits(),
                    narrow.swept(id)[j].to_bits(),
                    "{name}[{j}]"
                );
                assert!(
                    is_the_nearest_binary32(n[j], &row["series"][j])?,
                    "{name}[{j}]"
                );
                checked += 1;
            }
        }
        assert_eq!(checked, 8 * TERMS);
        Ok(())
    }

    #[test]
    fn the_helicoid_file_is_laid_out_as_rustfmt_leaves_it() -> Result<(), String> {
        let (wide, narrow) = series()?;
        let jsonl = series_file()?;
        for placeholder in [false, true] {
            let csv = match placeholder {
                true => synthetic([8; 8]),
                false => committed_csv()?,
            };
            let text = render_helicoid(&csv, &jsonl, (&wide, &narrow), placeholder)?;
            assert_eq!(rustfmt_diff(&text)?, None, "placeholder {placeholder}");
        }
        Ok(())
    }

    #[test]
    fn a_binary32_literal_is_rounded_from_the_rational_and_not_from_the_binary64_literal(
    ) -> Result<(), String> {
        // `k`'s first term made `(2^70 + 2^46 + 1)/2^70`, `1 + 2^-24 + 2^-70`: above the binary32
        // tie between 1 and `1 + 2^-23`, which binary64 holds as the tie itself and a second
        // rounding sends to the even 1. No term of the committed series is such a value, so the
        // committed file cannot show a literal taken from the binary64 one; this one does.
        let one = BigUint::from(1u32);
        let den = &one << 70usize;
        let term = format!("{}/{den}", &den + (&one << 46usize) + &one);
        let file = String::from_utf8(series_file()?).map_err(|e| e.to_string())?;
        let hazard = file.replacen(
            "\"series\":[\"1/2\",",
            &format!("\"series\":[\"{term}\","),
            1,
        );
        assert_ne!(hazard, file);
        let scratch = Scratch::new("hazard-series");
        std::fs::create_dir_all(&scratch.0).map_err(|e| e.to_string())?;
        std::fs::write(scratch.0.join(SERIES_FILE), &hazard).map_err(|e| e.to_string())?;
        let (wide, narrow) = (
            Series::<f64>::load(&scratch.0)?,
            Series::<f32>::load(&scratch.0)?,
        );
        let text = render(&committed_csv()?, hazard.as_bytes(), &wide, &narrow)?;
        let first = |name: &str| -> Result<String, String> {
            let literals = literals_of(&text, name)?;
            literals
                .first()
                .map(|l| (*l).to_string())
                .ok_or("no literal".into())
        };
        let (k32, k64) = (first("K_F32")?, first("K_F64")?);
        let (once, twice) = (
            k32.parse::<f32>().map_err(|e| e.to_string())?,
            k64.parse::<f64>().map_err(|e| e.to_string())? as f32,
        );
        assert_eq!(once.to_bits(), (1.0 + f32::EPSILON).to_bits(), "{k32}");
        assert_eq!(twice.to_bits(), 1.0f32.to_bits(), "{k64}");
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
            let terms = terms?;
            text.push_str(&around(&terms));
            // The same series as a constant, whose head is short or so long that a one-line array
            // within `array_width` still overflows `max_width`.
            for name in ["X", "SWEPT_COS_HALF_F64"] {
                let head = format!("pub(crate) const {name}: [f64; {}] = [", terms.len());
                text.push_str(&format!("{}\n", array_literal(&head, &terms, 0, "];")));
            }
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
        let (jsonl, (wide, narrow)) = (series_file()?, series()?);
        let base = synthetic([8; 8]);
        let refuses = |csv: &str, why: &str| {
            let e = render(csv, &jsonl, &wide, &narrow)
                .err()
                .unwrap_or_default();
            assert!(e.contains(why), "{why}: {e}");
        };
        let edit = |from: &str, to: &str| {
            assert!(base.contains(from), "{from}");
            base.replacen(from, to, 1)
        };
        refuses("", "documented header");
        refuses(&base.replacen("coeff,", "coeff ,", 1), "documented header");
        let first = base.lines().take(2).collect::<Vec<_>>().join("\n") + "\n";
        refuses(&first, "expected one per coefficient and precision, 16");
        refuses(&edit("c,f64,8", "c,f32,8"), "f64 row of `c`");
        refuses(&edit("k,f64,8", "c,f64,8"), "f64 row of `k`");
        refuses(&edit("r,f32,8", "r,f64,8"), "f32 row of `r`");
        refuses(&edit("c,f64,8", "c,f64,17"), "17 terms");
        refuses(&edit("c,f64,8", "c,f64,0"), "0 terms");
        refuses(&edit("c,f64,8", "c,f64,x"), "terms");
        // 2^-2 is 0x3fd0.. at `f64` and 0x3e80_0000 at `f32`; the decimal column must be the bits'
        // own, a z below 0 is refused, and so are the wrong number of digits for the precision.
        refuses(
            &edit("0x3fd0000000000000", "0x3fd0000000000001"),
            "disagree",
        );
        refuses(&edit("0x3e800000", "0x3e800001"), "disagree");
        refuses(
            &edit("0x3fd0000000000000", "3fd0000000000000"),
            "16 hex digits",
        );
        refuses(&edit("0x3e800000", "0x3fd0000000000000"), "8 hex digits");
        refuses(&edit("0x3fd0000000000000", "0x3e800000"), "16 hex digits");
        refuses(
            &edit("0x3fd0000000000000", "0xbfd0000000000000"),
            "disagree",
        );
        refuses(&edit("0x3e800000", "0xbe800000"), "disagree");
        refuses(&edit("1.5e0", "inf"), "not a finite number");
        refuses(&edit("2.5e1", "x"), "not a finite number");
        Ok(())
    }
}
