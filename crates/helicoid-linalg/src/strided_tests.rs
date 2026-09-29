//! Tests of `Strided` and `StridedMut`. Layouts are checked against the defining formula
//! `data[r * rs + c * cs]`, blocks against the same formula shifted, and every out-of-bounds
//! access must panic: a misfit view (which the constructors refuse in debug) is built through
//! the private fields, so the release path is tested in a debug build too. `overflow_never_wraps`
//! compares every access with an exact `u128` model over strides and indices up to `usize::MAX`.

use super::{fits, Strided, StridedMut};
use proptest::prelude::*;
use std::panic::{catch_unwind, AssertUnwindSafe};

/// `0, 1, 2, ...` as the backing store, so a cell's value is its index.
fn iota<const L: usize>() -> [u32; L] {
    core::array::from_fn(|i| i as u32)
}

/// Every `(r, c)` of `v` against `data[r * rs + c * cs]`.
fn check_layout(v: &Strided<'_, u32>, data: &[u32], rs: usize, cs: usize) {
    for r in 0..v.rows() {
        for c in 0..v.cols() {
            assert_eq!(v.get(r, c), data[r * rs + c * cs], "({r}, {c})");
        }
    }
}

#[test]
fn col_major_and_row_major_layouts() {
    let d = iota::<6>();
    let cm = Strided::col_major(&d, 3, 2);
    assert_eq!((cm.rows(), cm.cols()), (3, 2));
    check_layout(&cm, &d, 1, 3);
    assert_eq!((cm.get(2, 0), cm.get(0, 1)), (2, 3));
    let rm = Strided::row_major(&d, 3, 2);
    check_layout(&rm, &d, 2, 1);
    assert_eq!((rm.get(2, 0), rm.get(0, 1)), (4, 1));

    let mut m = iota::<6>();
    let mut w = StridedMut::col_major(&mut m, 3, 2);
    w.set(2, 1, 100);
    assert_eq!(
        (w.rows(), w.cols(), w.get(2, 1), w.get(0, 0)),
        (3, 2, 100, 0)
    );
    assert_eq!(m, [0, 1, 2, 3, 4, 100]);
    let mut m = iota::<6>();
    StridedMut::row_major(&mut m, 3, 2).set(2, 1, 100);
    assert_eq!(m, [0, 1, 2, 3, 4, 100]);
    let mut m = iota::<6>();
    StridedMut::row_major(&mut m, 3, 2).set(0, 1, 100);
    assert_eq!(m, [0, 100, 2, 3, 4, 5]);
}

#[test]
fn f64_entries_round_trip_bitwise() {
    let mut m = [0.0_f64; 4];
    let mut w = StridedMut::col_major(&mut m, 2, 2);
    for (r, c, v) in [
        (0, 0, -0.0),
        (1, 0, f64::NAN),
        (0, 1, 1e-310),
        (1, 1, f64::INFINITY),
    ] {
        w.set(r, c, v);
        assert_eq!(w.get(r, c).to_bits(), v.to_bits());
    }
    assert_eq!(
        Strided::col_major(&m, 2, 2).get(1, 0).to_bits(),
        f64::NAN.to_bits()
    );
}

/// The mapping of `with_strides` for a faer `MatMut` (column-major with leading dimension `ld`)
/// and a Ceres Jacobian buffer (row-major, `ld` columns): a 3x4 matrix in a padded buffer.
#[test]
fn faer_and_ceres_mappings() {
    let ld = 5;
    let d: [u32; 20] = iota();
    // faer: rs = 1, cs = ld; the last (ld - rows) entries of each column are padding.
    let faer = Strided::with_strides(&d, 3, 4, 1, ld);
    check_layout(&faer, &d, 1, ld);
    assert_eq!(faer.get(2, 3), 2 + 3 * 5);
    // Ceres: rs = ld, cs = 1; the last (ld - cols) entries of each row are padding.
    let ceres = Strided::with_strides(&d, 3, 4, ld, 1);
    check_layout(&ceres, &d, ld, 1);
    assert_eq!(ceres.get(2, 3), 2 * 5 + 3);
    // The transpose is the same data with the strides swapped.
    let t = Strided::with_strides(&d, 4, 3, ld, 1);
    for r in 0..3 {
        for c in 0..4 {
            assert_eq!(t.get(c, r), d[c * ld + r]);
            assert_eq!(t.get(c, r), faer.get(r, c));
        }
    }

    // A write through the padded view leaves the padding alone.
    let mut m = [0_u32; 20];
    let mut w = StridedMut::with_strides(&mut m, 3, 4, ld, 1);
    for r in 0..3 {
        for c in 0..4 {
            w.set(r, c, 1);
        }
    }
    let padding: u32 = (0..20).filter(|i| i % ld == 4).map(|i| m[i]).sum();
    assert_eq!((m.iter().sum::<u32>(), padding), (12, 0));
}

#[test]
fn blocks_and_blocks_of_blocks() {
    // A 4x5 row-major Jacobian; take rows 1..4, cols 1..5, then rows 1..3, cols 2..4 of that.
    let d: [u32; 20] = iota();
    let j = Strided::row_major(&d, 4, 5);
    let b = j.block(1, 1, 3, 4);
    let bb = b.block(1, 2, 2, 2);
    assert_eq!((bb.rows(), bb.cols()), (2, 2));
    for r in 0..2 {
        for c in 0..2 {
            assert_eq!(bb.get(r, c), j.get(r + 2, c + 3));
            assert_eq!(bb.get(r, c), d[(r + 2) * 5 + c + 3]);
        }
    }

    let mut m: [u32; 20] = iota();
    let mut w = StridedMut::row_major(&mut m, 4, 5);
    {
        let mut b = w.block(1, 1, 3, 4);
        let mut bb = b.block(1, 2, 2, 2);
        bb.set(1, 1, 900);
        assert_eq!(bb.get(1, 1), 900);
    }
    assert_eq!(w.get(3, 4), 900);
    assert_eq!(m[3 * 5 + 4], 900);
    assert_eq!(
        m.iter()
            .zip(iota::<20>())
            .filter(|(a, b)| **a != *b)
            .count(),
        1
    );
}

/// A write into a block changes exactly the block's cells; the sibling block is untouched.
#[test]
fn sibling_blocks_are_disjoint() {
    let mut m = [0_u32; 12];
    let mut w = StridedMut::col_major(&mut m, 3, 4);
    for (k, (c0, tag)) in [(0, 1_u32), (2, 2)].into_iter().enumerate() {
        let mut b = w.block(0, c0, 3, 2);
        for r in 0..3 {
            for c in 0..2 {
                b.set(r, c, tag * 10 + k as u32);
            }
        }
    }
    assert_eq!(m, [10, 10, 10, 10, 10, 10, 21, 21, 21, 21, 21, 21]);
}

#[test]
fn blocks_of_a_full_and_an_empty_extent() {
    let d: [u32; 6] = iota();
    let v = Strided::col_major(&d, 2, 3);
    let whole = v.block(0, 0, 2, 3);
    check_layout(&whole, &d, 1, 2);
    // Empty blocks are allowed anywhere up to and including the far edge.
    for (r0, c0, r, c) in [(2, 0, 0, 3), (0, 3, 2, 0), (2, 3, 0, 0), (1, 1, 0, 2)] {
        let e = v.block(r0, c0, r, c);
        assert_eq!((e.rows(), e.cols()), (r, c));
    }
}

#[test]
fn zero_size_views() {
    let empty: [u32; 0] = [];
    for (r, c) in [(0, 0), (0, 7), (7, 0)] {
        let v = Strided::col_major(&empty, r, c);
        assert_eq!((v.rows(), v.cols()), (r, c));
        assert!(catch_unwind(|| v.get(0, 0)).is_err());
        let mut e = [0_u32; 0];
        let mut w = StridedMut::row_major(&mut e, r, c);
        assert!(catch_unwind(AssertUnwindSafe(|| w.set(0, 0, 1))).is_err());
        let b = w.block(0, 0, r, c);
        assert_eq!((b.rows(), b.cols()), (r, c));
    }
    // A zero-size view of a non-empty slice, with strides that would not fit anything.
    let d = [1_u32, 2];
    let v = Strided::with_strides(&d, 0, 3, usize::MAX, usize::MAX);
    assert_eq!(v.block(0, 3, 0, 0).cols(), 0);
    // A 1x1 view needs only its origin, whatever the strides.
    let v = Strided::with_strides(&d, 1, 1, usize::MAX, usize::MAX);
    assert_eq!(v.get(0, 0), 1);
}

/// Strides of 0 broadcast in a read-only view; an overlapping mutable view is memory-safe and the
/// last write wins (documented on the type).
#[test]
fn broadcast_and_overlap() {
    let d = [7_u32, 8, 9];
    let row = Strided::with_strides(&d, 3, 3, 0, 1);
    assert_eq!((row.get(0, 2), row.get(2, 2), row.get(1, 0)), (9, 9, 7));

    let mut m = [0_u32; 3];
    let mut w = StridedMut::with_strides(&mut m, 2, 2, 1, 1);
    w.set(1, 0, 5);
    assert_eq!(w.get(0, 1), 5);
    w.set(0, 1, 6);
    assert_eq!(w.get(1, 0), 6);
    assert_eq!(m, [0, 6, 0]);
}

#[test]
fn fits_is_exact_and_overflow_safe() {
    assert!(fits(6, 3, 2, 1, 3) && !fits(5, 3, 2, 1, 3));
    assert!(fits(0, 0, 9, 1, 1) && fits(0, 9, 0, 1, 1) && !fits(0, 1, 1, 1, 1));
    assert!(fits(1, 1, 1, usize::MAX, usize::MAX));
    // (rows - 1) * rs overflows.
    assert!(!fits(usize::MAX, 3, 1, usize::MAX, 0));
    // Each product fits but their sum does not.
    assert!(!fits(usize::MAX, 2, 2, usize::MAX, usize::MAX));
    assert!(!fits(usize::MAX, usize::MAX, 2, 1, usize::MAX));
    // rows * cols overflows: a column-major usize::MAX x 2 view cannot fit.
    assert!(!fits(usize::MAX, usize::MAX, 2, 1, usize::MAX));
}

mod out_of_bounds {
    use super::*;

    #[test]
    #[should_panic(expected = "Strided::get: index out of range")]
    fn read_row_past_the_view() {
        let d = [0_u32; 6];
        // Row 3 of a 3x2 view is still inside the slice (index 3), but not in the view.
        let _ = Strided::col_major(&d, 3, 2).get(3, 0);
    }

    #[test]
    #[should_panic(expected = "Strided::get: index out of range")]
    fn read_column_past_the_view() {
        let d = [0_u32; 6];
        let _ = Strided::col_major(&d, 3, 2).get(0, 2);
    }

    #[test]
    #[should_panic(expected = "StridedMut::set: index out of range")]
    fn write_row_past_the_view() {
        let mut d = [0_u32; 6];
        StridedMut::col_major(&mut d, 3, 2).set(3, 0, 1);
    }

    #[test]
    #[should_panic(expected = "StridedMut::set: index out of range")]
    fn write_column_past_the_view() {
        let mut d = [0_u32; 6];
        StridedMut::row_major(&mut d, 3, 2).set(0, 2, 1);
    }

    /// A write past a block's last row would land in the next block of the same slice.
    #[test]
    #[should_panic(expected = "StridedMut::set: index out of range")]
    fn write_past_a_block_is_not_a_write_to_its_neighbour() {
        let mut d = [0_u32; 12];
        let mut w = StridedMut::col_major(&mut d, 3, 4);
        w.block(0, 0, 3, 2).set(0, 2, 1);
    }

    #[test]
    #[should_panic(expected = "StridedMut::get: index out of range")]
    fn mutable_read_row_past_the_view() {
        let mut d = [0_u32; 6];
        let _ = StridedMut::col_major(&mut d, 3, 2).get(3, 0);
    }

    #[test]
    #[should_panic(expected = "StridedMut::get: index out of range")]
    fn mutable_read_past_the_view() {
        let mut d = [0_u32; 6];
        let _ = StridedMut::col_major(&mut d, 3, 2).get(0, 2);
    }

    #[test]
    #[should_panic(expected = "Strided::block: block out of range")]
    fn block_past_the_rows() {
        let d = [0_u32; 6];
        let _ = Strided::col_major(&d, 3, 2).block(1, 0, 3, 1);
    }

    #[test]
    #[should_panic(expected = "StridedMut::block: block out of range")]
    fn mutable_block_past_the_columns() {
        let mut d = [0_u32; 6];
        let _ = StridedMut::col_major(&mut d, 3, 2).block(0, 1, 1, 2);
    }

    #[test]
    #[should_panic(expected = "StridedMut::block: block out of range")]
    fn empty_block_past_the_edge() {
        let mut d = [0_u32; 6];
        let _ = StridedMut::col_major(&mut d, 3, 2).block(4, 0, 0, 1);
    }

    #[test]
    #[should_panic(expected = "StridedMut::block: block out of range")]
    fn block_extent_overflow() {
        let mut d = [0_u32; 6];
        let _ = StridedMut::col_major(&mut d, 3, 2).block(2, 0, usize::MAX, 1);
    }

    /// A view whose strides reach past the slice, built past the constructor's debug check.
    fn misfit(d: &mut [u32]) -> StridedMut<'_, u32> {
        StridedMut {
            data: d,
            rows: 3,
            cols: 2,
            rs: 1,
            cs: 8,
        }
    }

    #[test]
    #[should_panic(expected = "out of")]
    fn write_past_the_slice() {
        let mut d = [0_u32; 6];
        misfit(&mut d).set(0, 1, 1);
    }

    #[test]
    #[should_panic(expected = "out of")]
    fn read_past_the_slice() {
        let mut d = [0_u32; 6];
        let _ = misfit(&mut d).get(2, 1);
    }

    #[test]
    #[should_panic(expected = "out of")]
    fn block_origin_past_the_slice() {
        let mut d = [0_u32; 6];
        let _ = misfit(&mut d).block(0, 1, 3, 1);
    }

    #[test]
    #[should_panic(expected = "StridedMut::block: origin out of the slice")]
    fn block_origin_at_the_slice_end() {
        // The first entry of column 1 is at index 6 == len: a legal empty tail, but no entry.
        let mut d = [0_u32; 6];
        let mut v = misfit(&mut d);
        v.cs = 6;
        let _ = v.block(0, 1, 3, 1);
    }

    #[test]
    #[should_panic(expected = "Strided::block: origin out of the slice")]
    fn read_block_origin_at_the_slice_end() {
        let d = [0_u32; 6];
        let v = Strided {
            data: &d,
            rows: 3,
            cols: 2,
            rs: 1,
            cs: 6,
        };
        let _ = v.block(0, 1, 3, 1);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "StridedMut: the view does not fit the slice")]
    fn constructor_checks_the_fit_in_debug() {
        let mut d = [0_u32; 5];
        let _ = StridedMut::col_major(&mut d, 3, 2);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "Strided: the view does not fit the slice")]
    fn read_constructor_checks_the_fit_in_debug() {
        let d = [0_u32; 6];
        let _ = Strided::with_strides(&d, 3, 2, 2, 3);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "does not fit")]
    fn constructor_check_survives_overflow() {
        let d = [0_u32; 6];
        let _ = Strided::col_major(&d, usize::MAX, 2);
    }
}

/// Widths and boundaries of `usize` at any pointer width; `2 * 2^(BITS - 1)` wraps to 0.
const BIG: [usize; 7] = [
    0,
    1,
    2,
    1 << (usize::BITS / 2),
    1 << (usize::BITS - 1),
    (1 << (usize::BITS - 1)) + 1,
    usize::MAX,
];

const LEN: usize = 4;
const CELLS: [u32; LEN] = [10, 11, 12, 13];

/// The exact index of `(r, c)` in a `LEN`-cell slice, `None` if outside the view or the slice.
fn exact(dims: (usize, usize), st: (usize, usize), r: usize, c: usize) -> Option<usize> {
    if r >= dims.0 || c >= dims.1 {
        return None;
    }
    (r as u128)
        .checked_mul(st.0 as u128)
        .zip((c as u128).checked_mul(st.1 as u128))
        .and_then(|(a, b)| a.checked_add(b))
        .filter(|&i| i < LEN as u128)
        .map(|i| i as usize)
}

/// `f`'s value, or `None` if it panicked. The panic must be an index or range failure of ours: a
/// debug-build arithmetic overflow ("attempt to multiply with overflow") is a defect, and only
/// the release profile would have hidden it.
fn attempt<T>(f: impl FnOnce() -> T) -> Option<T> {
    match catch_unwind(AssertUnwindSafe(f)) {
        Ok(x) => Some(x),
        Err(e) => {
            let msg = e
                .downcast_ref::<&str>()
                .copied()
                .or_else(|| e.downcast_ref::<std::string::String>().map(|s| s.as_str()))
                .unwrap_or("");
            assert!(
                msg.contains("out of") && !msg.contains("overflow"),
                "not an out-of-range panic: {msg:?}"
            );
            None
        }
    }
}

/// `f(rows, cols, rs, cs, i, j)` over every combination of `BIG`.
fn grid(mut f: impl FnMut((usize, usize), (usize, usize), usize, usize)) {
    for rows in BIG {
        for cols in BIG {
            for rs in BIG {
                for cs in BIG {
                    for i in BIG {
                        for j in BIG {
                            f((rows, cols), (rs, cs), i, j);
                        }
                    }
                }
            }
        }
    }
}

fn mut_view(buf: &mut [u32; LEN], d: (usize, usize), st: (usize, usize)) -> StridedMut<'_, u32> {
    StridedMut {
        data: buf,
        rows: d.0,
        cols: d.1,
        rs: st.0,
        cs: st.1,
    }
}

fn read_view(buf: &[u32; LEN], d: (usize, usize), st: (usize, usize)) -> Strided<'_, u32> {
    Strided {
        data: buf,
        rows: d.0,
        cols: d.1,
        rs: st.0,
        cs: st.1,
    }
}

/// Every access of a misfit view, through both view types, either panics or agrees with the
/// exact `u128` index; a stride product that wraps modulo 2^BITS to a small index panics instead.
/// `get` runs alone, so its own row and column checks are what refuse an out-of-view read.
#[test]
fn overflow_never_wraps() {
    let (mut in_bounds, mut refused) = (0_u32, 0_u32);
    grid(|d, st, r, c| {
        let want = exact(d, st, r, c);

        let mut buf = CELLS;
        let got = attempt(|| mut_view(&mut buf, d, st).get(r, c));
        assert_eq!(
            got,
            want.map(|i| CELLS[i]),
            "mut get {d:?} {st:?} ({r}, {c})"
        );
        let got = attempt(|| read_view(&CELLS, d, st).get(r, c));
        assert_eq!(got, want.map(|i| CELLS[i]), "get {d:?} {st:?} ({r}, {c})");

        let mut buf = CELLS;
        let mut v = mut_view(&mut buf, d, st);
        let done = attempt(|| {
            v.set(r, c, 99);
            v.get(r, c)
        });
        let mut expect = CELLS;
        if let Some(i) = want {
            expect[i] = 99;
            in_bounds += 1;
        } else {
            refused += 1;
        }
        assert_eq!(done, want.map(|_| 99), "set {d:?} {st:?} ({r}, {c})");
        assert_eq!(buf, expect, "set {d:?} {st:?} ({r}, {c})");
    });
    assert!(in_bounds > 1000 && refused > 1000, "{in_bounds} {refused}");
}

/// The block path on the same grid, for both view types: a block is refused, or its origin is the
/// exact index (an origin at `LEN` included, which is a legal empty tail but no entry).
#[test]
fn block_overflow_never_wraps() {
    grid(|d, st, r0, c0| {
        let want = exact(d, st, r0, c0).map(|i| CELLS[i]);
        let mut buf = CELLS;
        let got = attempt(|| mut_view(&mut buf, d, st).block(r0, c0, 1, 1).get(0, 0));
        assert_eq!(got, want, "mut block {d:?} {st:?} ({r0}, {c0})");
        let got = attempt(|| read_view(&CELLS, d, st).block(r0, c0, 1, 1).get(0, 0));
        assert_eq!(got, want, "block {d:?} {st:?} ({r0}, {c0})");
    });
}

proptest! {
    #![proptest_config(ProptestConfig { failure_persistence: None, ..ProptestConfig::default() })]

    /// A view that fits, a block of it and a block of that block read `data[r * rs + c * cs]`
    /// shifted by the block origins; a write changes that one cell and no other, unless two cells
    /// share an index, when the write is visible at both.
    #[test]
    fn views_match_the_index_formula(
        rows in 0_usize..6, cols in 0_usize..6, rs in 0_usize..8, cs in 0_usize..8,
        pad in 0_usize..3, seed in any::<u32>(),
        sel in proptest::array::uniform8(any::<usize>()),
    ) {
        let need = if rows == 0 || cols == 0 { 0 } else { (rows - 1) * rs + (cols - 1) * cs + 1 };
        let len = need + pad;
        prop_assert!(fits(len, rows, cols, rs, cs));
        let d: std::vec::Vec<u32> = (0..len as u32).map(|i| i ^ seed).collect();
        let v = Strided::with_strides(&d, rows, cols, rs, cs);
        // A block `(r0, c0, r, c)` of an `h x w` view, from four selectors.
        let pick = |s: &[usize], h: usize, w: usize| {
            let (r0, c0) = (s[0] % (h + 1), s[1] % (w + 1));
            (r0, c0, s[2] % (h - r0 + 1), s[3] % (w - c0 + 1))
        };
        let b = pick(&sel[..4], rows, cols);
        let blk = v.block(b.0, b.1, b.2, b.3);
        let bb = pick(&sel[4..], b.2, b.3);
        let blk2 = blk.block(bb.0, bb.1, bb.2, bb.3);
        for r in 0..bb.2 {
            for c in 0..bb.3 {
                let (rr, cc) = (b.0 + bb.0 + r, b.1 + bb.1 + c);
                prop_assert_eq!(blk2.get(r, c), d[rr * rs + cc * cs]);
                prop_assert_eq!(blk2.get(r, c), v.get(rr, cc));
            }
        }

        let mut m = d.clone();
        let mut w = StridedMut::with_strides(&mut m, rows, cols, rs, cs);
        let mut wb = w.block(b.0, b.1, b.2, b.3);
        let mut wb2 = wb.block(bb.0, bb.1, bb.2, bb.3);
        for r in 0..bb.2 {
            for c in 0..bb.3 {
                wb2.set(r, c, u32::MAX);
            }
        }
        let touched: std::collections::BTreeSet<usize> = (0..bb.2)
            .flat_map(|r| (0..bb.3).map(move |c| (b.0 + bb.0 + r) * rs + (b.1 + bb.1 + c) * cs))
            .collect();
        for (i, (&after, &before)) in m.iter().zip(&d).enumerate() {
            prop_assert_eq!(after, if touched.contains(&i) { u32::MAX } else { before });
        }
    }
}
