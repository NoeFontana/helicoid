//! Strided views over caller memory (`docs/PHASE2.md` §5).
//!
//! A view maps `(r, c)` to `data[r * rs + c * cs]`. Solvers hand `helicoid` a block of their
//! Jacobian through it: the usual column-major faer `MatMut` is
//! `with_strides(data, rows, cols, 1, col_stride)` and a Ceres row-major buffer is `row_major`, or
//! `with_strides(data, rows, cols, ld, 1)` for a padded one. Strides are non-negative: a reversed
//! (negative-stride) view is copied or re-based by the caller.
//!
//! Every access checks `r < rows`, `c < cols` and the index against the slice, in release too:
//! an out-of-bounds access panics (D11's one documented panic class), and so does a `block` whose
//! first entry is past the slice. The index is computed with
//! saturating arithmetic, so a product that overflows `usize` becomes `usize::MAX`, which is
//! past the end of every slice, and never wraps to an in-bounds cell. The row and column checks
//! keep a write inside its own view: without them a write past the last row of a block would
//! land in its neighbour.

/// `true` iff every `(r, c)` with `r < rows`, `c < cols` maps below `len`, without overflow.
fn fits(len: usize, rows: usize, cols: usize, rs: usize, cs: usize) -> bool {
    if rows == 0 || cols == 0 {
        return true;
    }
    let last = (rows - 1)
        .checked_mul(rs)
        .zip((cols - 1).checked_mul(cs))
        .and_then(|(a, b)| a.checked_add(b));
    last.is_some_and(|i| i < len)
}

/// The slice index of `(r, c)`, `usize::MAX` if the arithmetic overflows.
#[inline]
fn index(rs: usize, cs: usize, r: usize, c: usize) -> usize {
    r.saturating_mul(rs).saturating_add(c.saturating_mul(cs))
}

/// `true` iff the block `(r0, c0, rows, cols)` lies inside a `have_rows x have_cols` view.
#[inline]
fn block_fits(have: (usize, usize), r0: usize, c0: usize, rows: usize, cols: usize) -> bool {
    r0.checked_add(rows).is_some_and(|e| e <= have.0)
        && c0.checked_add(cols).is_some_and(|e| e <= have.1)
}

/// A read-only `rows x cols` view of a slice with row stride `rs` and column stride `cs`.
///
/// `Copy`: a view is a borrow. Layout is not a contract (D2).
#[derive(Clone, Copy, Debug)]
pub struct Strided<'a, S> {
    data: &'a [S],
    rows: usize,
    cols: usize,
    rs: usize,
    cs: usize,
}

impl<'a, S: Copy> Strided<'a, S> {
    /// A column-major `rows x cols` view: entry `(r, c)` is `data[r + c * rows]`.
    ///
    /// # Domain
    ///
    /// The view fits the slice, without overflow, checked by `debug_assert!`; see
    /// [`Strided::with_strides`].
    #[inline]
    pub fn col_major(data: &'a [S], rows: usize, cols: usize) -> Self {
        Self::with_strides(data, rows, cols, 1, rows)
    }

    /// A row-major `rows x cols` view: entry `(r, c)` is `data[r * cols + c]`.
    ///
    /// # Domain
    ///
    /// The view fits the slice, without overflow, checked by `debug_assert!`; see
    /// [`Strided::with_strides`].
    #[inline]
    pub fn row_major(data: &'a [S], rows: usize, cols: usize) -> Self {
        Self::with_strides(data, rows, cols, cols, 1)
    }

    /// A view with entry `(r, c)` at `data[r * rs + c * cs]`. Strides are in elements. A stride
    /// of 0 broadcasts, and overlapping strides are allowed for a read-only view.
    ///
    /// # Domain
    ///
    /// The last entry fits: for `rows, cols > 0`, `(rows - 1) * rs + (cols - 1) * cs < data.len()`
    /// without overflow, checked by `debug_assert!`. A release build does not check here; an
    /// access to an entry outside the slice panics (see [`Strided::get`]).
    #[inline]
    pub fn with_strides(data: &'a [S], rows: usize, cols: usize, rs: usize, cs: usize) -> Self {
        debug_assert!(
            fits(data.len(), rows, cols, rs, cs),
            "Strided: the view does not fit the slice"
        );
        Self {
            data,
            rows,
            cols,
            rs,
            cs,
        }
    }

    /// The number of rows.
    #[inline]
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// The number of columns.
    #[inline]
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Entry `(r, c)`.
    ///
    /// # Panics
    ///
    /// If `r >= rows`, `c >= cols` or the entry lies outside the slice, in release too.
    #[inline]
    pub fn get(&self, r: usize, c: usize) -> S {
        assert!(
            r < self.rows && c < self.cols,
            "Strided::get: index out of range"
        );
        self.data[index(self.rs, self.cs, r, c)]
    }

    /// The `rows x cols` block whose top-left entry is `(r0, c0)`; same strides, same slice.
    ///
    /// # Panics
    ///
    /// If `r0 + rows > self.rows()` or `c0 + cols > self.cols()` (also for an empty block), or the
    /// block is non-empty and its first entry lies outside the slice.
    #[inline]
    pub fn block(&self, r0: usize, c0: usize, rows: usize, cols: usize) -> Strided<'a, S> {
        assert!(
            block_fits((self.rows, self.cols), r0, c0, rows, cols),
            "Strided::block: block out of range"
        );
        let data: &'a [S] = self.data;
        let data = if rows == 0 || cols == 0 {
            &data[..0]
        } else {
            let i = index(self.rs, self.cs, r0, c0);
            assert!(i < data.len(), "Strided::block: origin out of the slice");
            &data[i..]
        };
        Strided {
            data,
            rows,
            cols,
            rs: self.rs,
            cs: self.cs,
        }
    }
}

/// A writable `rows x cols` view of a slice with row stride `rs` and column stride `cs`.
///
/// Not `Copy`: it holds the unique borrow. Layout is not a contract (D2).
///
/// Overlapping strides (two entries at one index, e.g. `rs = cs = 1` on a `2 x 2` view) are not
/// rejected: they are memory-safe and the last write wins, but the view is then no longer a
/// matrix. Solver layouts (dense, padded, transposed) never overlap.
#[derive(Debug)]
pub struct StridedMut<'a, S> {
    data: &'a mut [S],
    rows: usize,
    cols: usize,
    rs: usize,
    cs: usize,
}

impl<'a, S: Copy> StridedMut<'a, S> {
    /// A column-major `rows x cols` view: entry `(r, c)` is `data[r + c * rows]`.
    ///
    /// # Domain
    ///
    /// The view fits the slice, without overflow, checked by `debug_assert!`; see
    /// [`StridedMut::with_strides`].
    #[inline]
    pub fn col_major(data: &'a mut [S], rows: usize, cols: usize) -> Self {
        Self::with_strides(data, rows, cols, 1, rows)
    }

    /// A row-major `rows x cols` view: entry `(r, c)` is `data[r * cols + c]`.
    ///
    /// # Domain
    ///
    /// The view fits the slice, without overflow, checked by `debug_assert!`; see
    /// [`StridedMut::with_strides`].
    #[inline]
    pub fn row_major(data: &'a mut [S], rows: usize, cols: usize) -> Self {
        Self::with_strides(data, rows, cols, cols, 1)
    }

    /// A view with entry `(r, c)` at `data[r * rs + c * cs]`. Strides are in elements.
    ///
    /// # Domain
    ///
    /// The last entry fits: for `rows, cols > 0`, `(rows - 1) * rs + (cols - 1) * cs < data.len()`
    /// without overflow, checked by `debug_assert!`. A release build does not check here; an
    /// access to an entry outside the slice panics (see [`StridedMut::set`]). Strides that make
    /// two entries share an index are not checked (see the type's documentation).
    #[inline]
    pub fn with_strides(data: &'a mut [S], rows: usize, cols: usize, rs: usize, cs: usize) -> Self {
        debug_assert!(
            fits(data.len(), rows, cols, rs, cs),
            "StridedMut: the view does not fit the slice"
        );
        Self {
            data,
            rows,
            cols,
            rs,
            cs,
        }
    }

    /// The number of rows.
    #[inline]
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// The number of columns.
    #[inline]
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Entry `(r, c)`.
    ///
    /// # Panics
    ///
    /// If `r >= rows`, `c >= cols` or the entry lies outside the slice, in release too.
    #[inline]
    pub fn get(&self, r: usize, c: usize) -> S {
        assert!(
            r < self.rows && c < self.cols,
            "StridedMut::get: index out of range"
        );
        self.data[index(self.rs, self.cs, r, c)]
    }

    /// Set entry `(r, c)` to `v`.
    ///
    /// # Panics
    ///
    /// If `r >= rows`, `c >= cols` or the entry lies outside the slice, in release too: an
    /// out-of-bounds write never wraps around and never lands outside the view (D11).
    #[inline]
    pub fn set(&mut self, r: usize, c: usize, v: S) {
        assert!(
            r < self.rows && c < self.cols,
            "StridedMut::set: index out of range"
        );
        self.data[index(self.rs, self.cs, r, c)] = v;
    }

    /// The `rows x cols` block whose top-left entry is `(r0, c0)`; same strides, same slice. It
    /// borrows `self` mutably, and a block of a block composes.
    ///
    /// # Panics
    ///
    /// If `r0 + rows > self.rows()` or `c0 + cols > self.cols()` (also for an empty block), or the
    /// block is non-empty and its first entry lies outside the slice.
    #[inline]
    pub fn block(&mut self, r0: usize, c0: usize, rows: usize, cols: usize) -> StridedMut<'_, S> {
        assert!(
            block_fits((self.rows, self.cols), r0, c0, rows, cols),
            "StridedMut::block: block out of range"
        );
        let data = if rows == 0 || cols == 0 {
            &mut self.data[..0]
        } else {
            let i = index(self.rs, self.cs, r0, c0);
            assert!(
                i < self.data.len(),
                "StridedMut::block: origin out of the slice"
            );
            &mut self.data[i..]
        };
        StridedMut {
            data,
            rows,
            cols,
            rs: self.rs,
            cs: self.cs,
        }
    }
}

#[cfg(test)]
#[path = "strided_tests.rs"]
mod tests;
