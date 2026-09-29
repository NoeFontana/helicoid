//! The search of `docs/PHASE1.md` §6 over per-record errors formed once. A candidate is a number
//! of series terms and a switch; on a record it uses the series arm below the switch (`z <
//! switch`, the kernel's own comparison) and the exact arm at or above it, so its score is a
//! maximum over records of one arm's error each: exact, and cheap for every candidate.

use std::cmp::Ordering;

/// Series lengths swept, `1..=TERMS`.
pub(super) const TERMS: usize = 8;

/// The error of one arm at one record, in units of `u`: the value and `d/dz`. Infinite when the
/// output is not finite (`metric::Score::NonFinite`), so an arm that cannot be selected safely
/// makes every candidate that selects it infinite.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Errors {
    pub(super) value: f64,
    pub(super) deriv: f64,
}

/// One record: its branch variable and the error of each arm there.
#[derive(Clone, Copy, Debug)]
pub(super) struct Sample {
    pub(super) z: f64,
    pub(super) exact: Errors,
    /// `series[m - 1]` is the arm of `m` terms.
    pub(super) series: [Errors; TERMS],
}

/// The field of a record's error that a maximum is over.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Field {
    Value,
    Deriv,
}

impl Field {
    /// The name in the CSV, as in `value_max_u`.
    pub(super) fn name(self) -> &'static str {
        match self {
            Field::Value => "value",
            Field::Deriv => "deriv",
        }
    }
}

/// A candidate's maxima over the records, and where each is attained (an index into the
/// samples, the first on a tie).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Score {
    pub(super) terms: usize,
    pub(super) switch: f64,
    pub(super) value: f64,
    pub(super) deriv: f64,
    value_at: usize,
    deriv_at: usize,
}

impl Score {
    /// The objective of `0004` item 1: the larger of the two maxima.
    pub(super) fn objective(&self) -> f64 {
        self.value.max(self.deriv)
    }

    /// The field and the record that attain the objective: the value on a tie.
    pub(super) fn argmax(&self) -> (Field, usize) {
        match self.value >= self.deriv {
            true => (Field::Value, self.value_at),
            false => (Field::Deriv, self.deriv_at),
        }
    }
}

/// `terms` series terms below `switch`, over `samples`.
pub(super) fn score(samples: &[Sample], terms: usize, switch: f64) -> Score {
    debug_assert!((1..=TERMS).contains(&terms));
    let mut score = Score {
        terms,
        switch,
        value: 0.0,
        deriv: 0.0,
        value_at: 0,
        deriv_at: 0,
    };
    for (i, s) in samples.iter().enumerate() {
        let e = if s.z < switch {
            s.series[terms - 1]
        } else {
            s.exact
        };
        if e.value > score.value {
            (score.value, score.value_at) = (e.value, i);
        }
        if e.deriv > score.deriv {
            (score.deriv, score.deriv_at) = (e.deriv, i);
        }
    }
    score
}

/// Lowest objective first; ties to fewer terms, then to the larger switch (`docs/PHASE1.md` §6).
fn preference(a: &Score, b: &Score) -> Ordering {
    let by = a.objective().total_cmp(&b.objective());
    by.then(a.terms.cmp(&b.terms))
        .then(b.switch.total_cmp(&a.switch))
}

/// The outcome of one search.
#[derive(Debug)]
pub(super) struct Sweep {
    pub(super) chosen: Score,
    /// The position of the chosen switch in the grid.
    pub(super) index: usize,
    /// The best objective among the candidates at the top of the grid.
    pub(super) top: f64,
    pub(super) prior: Score,
    /// 1 + the grid candidates with a strictly smaller objective than the prior's.
    pub(super) prior_rank: usize,
    /// Grid candidates with exactly the chosen objective: the tie-break chose among them.
    pub(super) tied: usize,
    /// The smallest grid objective above the chosen one; `None` when every candidate ties.
    pub(super) next: Option<f64>,
    /// The chosen terms one grid point below and above the chosen switch; `None` at an end.
    pub(super) below: Option<f64>,
    pub(super) above: Option<f64>,
}

/// Every `(terms, switch)` of `1..=TERMS` times `grid` (increasing), and the prior candidate
/// `(terms, switch)`, over `samples`.
pub(super) fn search(
    samples: &[Sample],
    grid: &[f64],
    prior: (usize, f64),
) -> Result<Sweep, String> {
    let mut all = Vec::with_capacity(TERMS * grid.len());
    for terms in 1..=TERMS {
        for (i, &switch) in grid.iter().enumerate() {
            all.push((i, score(samples, terms, switch)));
        }
    }
    let (index, chosen) = *all
        .iter()
        .min_by(|a, b| preference(&a.1, &b.1))
        .ok_or("the switch grid is empty")?;
    let prior = score(samples, prior.0, prior.1);
    let objectives = || all.iter().map(|(_, s)| s.objective());
    let (best, worse) = (chosen.objective(), prior.objective());
    let neighbour = |i: Option<usize>| {
        i.and_then(|i| grid.get(i))
            .map(|&switch| score(samples, chosen.terms, switch).objective())
    };
    Ok(Sweep {
        chosen,
        index,
        top: all
            .iter()
            .filter(|(i, _)| *i + 1 == grid.len())
            .map(|(_, s)| s.objective())
            .min_by(f64::total_cmp)
            .unwrap_or(f64::INFINITY),
        prior,
        prior_rank: 1 + objectives().filter(|o| o.total_cmp(&worse).is_lt()).count(),
        tied: objectives().filter(|o| o.total_cmp(&best).is_eq()).count(),
        next: objectives()
            .filter(|o| o.total_cmp(&best).is_gt())
            .min_by(f64::total_cmp),
        below: neighbour(index.checked_sub(1)),
        above: neighbour(Some(index + 1)),
    })
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;
    use crate::thresholds::grid::grid;

    fn at(value: f64, deriv: f64) -> Errors {
        Errors { value, deriv }
    }

    /// An error of `x` in the derivative and half of it in the value: the derivative decides.
    fn deriv(x: f64) -> Errors {
        at(x / 2.0, x)
    }

    /// A record at each grid point, with the errors `exact(z)` and `series(m, z)`.
    fn synthetic(
        grid: &[f64],
        exact: impl Fn(f64) -> Errors,
        series: impl Fn(usize, f64) -> Errors,
    ) -> Vec<Sample> {
        let sample = |&z: &f64| Sample {
            z,
            exact: exact(z),
            series: std::array::from_fn(|i| series(i + 1, z)),
        };
        grid.iter().map(sample).collect()
    }

    #[test]
    fn the_search_finds_a_known_optimum() -> Result<(), String> {
        let g = grid();
        // The exact arm is bad below `g[600]` and the series above it, so only a switch exactly
        // there survives; a term fewer than 8 costs one unit more.
        let t = g[600];
        let samples = synthetic(
            &g,
            |z| deriv(if z < t { 1e6 } else { 1.0 }),
            |m, z| deriv(if z < t { (9 - m) as f64 } else { 1e9 }),
        );
        let s = search(&samples, &g, (4, g[10]))?;
        assert_eq!((s.chosen.terms, s.index), (8, 600));
        assert_eq!(s.chosen.switch.to_bits(), t.to_bits());
        // The derivative is the larger field, and the maxima are each field's own.
        assert_eq!((s.chosen.value, s.chosen.deriv), (0.5, 1.0));
        assert_eq!((s.chosen.objective(), s.tied, s.next), (1.0, 1, Some(2.0)));
        // A switch one point either side puts one record on the wrong arm.
        assert_eq!((s.below, s.above), (Some(1e6), Some(1e9)));
        // The top of the grid keeps the series arm on every record but the last: 1e9.
        assert_eq!(s.top, 1e9);
        // The prior at `g[10]` selects the exact arm below `t`: 1e6, beaten by every candidate
        // that has the switch right, which is one per length.
        assert_eq!((s.prior.objective(), s.prior_rank), (1e6, 9));
        Ok(())
    }

    #[test]
    fn ties_go_to_fewer_terms_then_to_the_larger_switch() -> Result<(), String> {
        let g = grid();
        // The exact arm is right from `g[600]`, so a switch below it is bad. A series arm is
        // right up to `g[650]` for 3 terms and up to `g[700]` for 4 or more (1 and 2 terms cost 5).
        // The tied candidates are 3 terms with `g[600]..=g[651]` and 4 to 8 with `g[600]..=g[701]`:
        // the largest switch of all is 4 terms at 701, so only the terms decide first.
        let reach = |m| if m == 3 { g[650] } else { g[700] };
        let samples = synthetic(
            &g,
            |z| deriv(if z < g[600] { 1e6 } else { 1.0 }),
            |m, z| {
                deriv(match (z <= reach(m), m <= 2) {
                    (false, _) => 1e9,
                    (true, true) => 5.0,
                    (true, false) => 1.0,
                })
            },
        );
        let s = search(&samples, &g, (2, g[650]))?;
        assert_eq!((s.chosen.terms, s.index), (3, 651));
        let ties = 52 + 5 * 102;
        assert_eq!((s.tied, s.next), (ties, Some(5.0)));
        // The prior scores 5: beaten by every one of the tied candidates.
        assert_eq!((s.prior.objective(), s.prior_rank), (5.0, 1 + ties));
        Ok(())
    }

    #[test]
    fn everything_tied_has_no_next_objective_and_an_empty_grid_is_an_error() -> Result<(), String> {
        let g = grid();
        let flat = synthetic(&g, |_| deriv(3.0), |_, _| deriv(3.0));
        let s = search(&flat, &g, (4, 0.01))?;
        // One term and the largest switch: the preference alone decides.
        assert_eq!(
            (s.chosen.terms, s.index, s.tied, s.next),
            (1, 1024, 8 * 1025, None)
        );
        assert_eq!((s.below, s.above), (Some(3.0), None));
        assert_eq!(s.top, 3.0);
        assert!(search(&flat, &[], (4, 0.01)).is_err());
        Ok(())
    }

    #[test]
    fn a_candidate_selects_by_the_kernels_own_comparison() {
        let (bad, good) = (at(9.0, 0.0), at(1.0, 0.0));
        let sample = |z, exact, series| Sample {
            z,
            exact,
            series: [series; TERMS],
        };
        let two = [sample(0.5, bad, good), sample(1.0, good, bad)];
        // `z < switch` is strict: the record at the switch is the exact arm's.
        assert_eq!(score(&two, 3, 1.0).objective(), 1.0);
        assert_eq!(score(&two, 3, 1.5).objective(), 9.0);
        assert_eq!(score(&two, 3, 0.5).objective(), 9.0);
    }

    #[test]
    fn the_objective_is_the_larger_maximum_and_names_where_it_is_attained() {
        let sample = |z, exact, series| Sample {
            z,
            exact,
            series: [series; TERMS],
        };
        let s = [
            sample(0.1, at(1.0, 2.0), at(9.0, 3.0)),
            sample(0.5, at(6.0, 9.0), at(5.0, 8.0)),
            sample(2.0, at(3.0, 5.0), at(0.0, 0.0)),
            sample(3.0, at(3.0, 5.0), at(0.0, 0.0)),
        ];
        // Each maximum is over its own field; a tie between them names the value.
        let seen = |switch| {
            let c = score(&s, 2, switch);
            (c.value, c.deriv, c.objective(), c.argmax())
        };
        assert_eq!(seen(1.0), (9.0, 8.0, 9.0, (Field::Value, 0)));
        assert_eq!(seen(0.3), (9.0, 9.0, 9.0, (Field::Value, 0)));
        assert_eq!(seen(0.0), (6.0, 9.0, 9.0, (Field::Deriv, 1)));
    }

    #[test]
    fn an_infinite_error_never_wins_against_a_finite_one() -> Result<(), String> {
        let g = grid();
        let s = search(
            &synthetic(
                &g,
                |z| {
                    if z < g[5] {
                        at(1.0, f64::INFINITY)
                    } else {
                        deriv(2.0)
                    }
                },
                |_, _| deriv(2.0),
            ),
            &g,
            (4, g[5]),
        )?;
        assert_eq!((s.chosen.objective(), s.next), (2.0, Some(f64::INFINITY)));
        Ok(())
    }
}
