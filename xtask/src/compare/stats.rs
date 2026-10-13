//! The statistics of a paired comparison (`0038` decisions 2, 5 and 6; the readings 0038 leaves
//! open are `0078`'s): per record, which subject's error is smaller and by what ratio; per
//! group of records, the counts, an exact sign test, the geometric mean of the ratio and its two
//! extremes. Pure, so it tests without the corpus.

use std::cmp::Ordering;

use crate::conformance::metric::Score;

/// The paired tally of a group of records, `b` read against `a`: "better" means `b`'s error is
/// smaller. A group is a stratum, an id (its strata pooled) or the whole comparison (its ids
/// pooled), by [`Tally::merge`].
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct Tally {
    /// Records evaluated, every kind below included.
    pub(crate) n: usize,
    /// `b < a`.
    pub(crate) better: usize,
    /// `b > a`.
    pub(crate) worse: usize,
    /// `b == a`, to the bit of the error: out of the sign test and the ratios.
    pub(crate) ties: usize,
    /// Of `better + worse`, the records where one error is exactly zero: in the sign test, out of
    /// the geometric mean, whose logarithm is infinite there. The extremes read them as `inf`.
    pub(crate) zero_sided: usize,
    /// A non-finite output (or an overflowing error) on `a`'s side, `b`'s, or both: out of every
    /// statistic, since a non-finite output is a conformance failure and not a magnitude.
    pub(crate) nonfinite_a: usize,
    pub(crate) nonfinite_b: usize,
    /// `Σ ln(b/a)` and its count, over the records that differ with both errors positive.
    log_sum: f64,
    logs: usize,
    /// The largest `b/a` over the records where `b` is worse, and the largest `a/b` where it is
    /// better; `None` when there is no such record.
    pub(crate) worst_b: Option<f64>,
    pub(crate) worst_a: Option<f64>,
    /// The per-stratum maxima over the records both scored finite, as the bar reads each subject.
    /// Meaningful for a stratum; for a pooled group it is the maximum of the maxima.
    pub(crate) max_a: Option<f64>,
    pub(crate) max_b: Option<f64>,
    /// How many of the pooled strata have `max_b / max_a` below, above and at 1 (the hand tally of
    /// `0076` decision 1); zero for a single stratum's tally until [`Tally::close_stratum`].
    pub(crate) strata_better: usize,
    pub(crate) strata_worse: usize,
    pub(crate) strata_tied: usize,
    /// The smallest and the largest `max_b / max_a` over the pooled strata.
    pub(crate) strata_ratio: Option<(f64, f64)>,
}

fn max_opt(a: Option<f64>, b: f64) -> Option<f64> {
    Some(a.map_or(b, |m| m.max(b)))
}

impl Tally {
    /// One record's two scores.
    pub(crate) fn add(&mut self, a: Score, b: Score) {
        self.n += 1;
        // A `Finite` NaN is not a score the metric emits; it is read as non-finite, not compared.
        let finite = |s: Score| match s {
            Score::Finite(x) if x.is_finite() => Some(x),
            Score::Finite(_) | Score::NonFinite => None,
        };
        let (a, b) = match (finite(a), finite(b)) {
            (Some(a), Some(b)) => (a, b),
            (a, b) => {
                self.nonfinite_a += usize::from(a.is_none());
                self.nonfinite_b += usize::from(b.is_none());
                return;
            }
        };
        self.max_a = max_opt(self.max_a, a);
        self.max_b = max_opt(self.max_b, b);
        // An error is never negative, so `total_cmp`'s order is the value order: equal is a tie to
        // the bit, and one zero of two unequal errors is the smaller.
        match b.total_cmp(&a) {
            Ordering::Equal => {
                self.ties += 1;
                return;
            }
            Ordering::Less => {
                self.better += 1;
                self.worst_a = max_opt(self.worst_a, a / b);
            }
            Ordering::Greater => {
                self.worse += 1;
                self.worst_b = max_opt(self.worst_b, b / a);
            }
        }
        if a.min(b) <= 0.0 {
            self.zero_sided += 1;
        } else {
            self.log_sum += (b / a).ln();
            self.logs += 1;
        }
    }

    /// `max_b / max_a`: the bar's two numbers as one ratio. `1` when both are zero, `inf` when
    /// only `a`'s is; `None` when no record was scored on both sides.
    pub(crate) fn max_ratio(&self) -> Option<f64> {
        let (a, b) = (self.max_a?, self.max_b?);
        Some(match b.total_cmp(&a) {
            Ordering::Equal => 1.0,
            _ => b / a,
        })
    }

    /// Counts this tally, read as one stratum, into the per-stratum maximum tallies.
    pub(crate) fn close_stratum(&mut self) {
        self.strata_ratio = self.max_ratio().map(|r| (r, r));
        // By the maxima themselves, not their rounded quotient.
        if let (Some(a), Some(b)) = (self.max_a, self.max_b) {
            match b.total_cmp(&a) {
                Ordering::Less => self.strata_better += 1,
                Ordering::Greater => self.strata_worse += 1,
                Ordering::Equal => self.strata_tied += 1,
            }
        }
    }

    /// Pools `other` into `self`: counts add, sums add, extremes take the larger.
    pub(crate) fn merge(&mut self, other: &Tally) {
        self.n += other.n;
        self.better += other.better;
        self.worse += other.worse;
        self.ties += other.ties;
        self.zero_sided += other.zero_sided;
        self.nonfinite_a += other.nonfinite_a;
        self.nonfinite_b += other.nonfinite_b;
        self.log_sum += other.log_sum;
        self.logs += other.logs;
        let pick = |x: Option<f64>, y: Option<f64>| match (x, y) {
            (Some(x), Some(y)) => Some(x.max(y)),
            (x, y) => x.or(y),
        };
        self.worst_a = pick(self.worst_a, other.worst_a);
        self.worst_b = pick(self.worst_b, other.worst_b);
        self.max_a = pick(self.max_a, other.max_a);
        self.max_b = pick(self.max_b, other.max_b);
        self.strata_better += other.strata_better;
        self.strata_worse += other.strata_worse;
        self.strata_tied += other.strata_tied;
        self.strata_ratio = match (self.strata_ratio, other.strata_ratio) {
            (Some((lo, hi)), Some((l, h))) => Some((lo.min(l), hi.max(h))),
            (x, y) => x.or(y),
        };
    }

    /// The geometric mean of `b/a` over the records that differ with both errors positive;
    /// `None` when there is none.
    pub(crate) fn geomean(&self) -> Option<f64> {
        (self.logs > 0).then(|| (self.log_sum / self.logs as f64).exp())
    }

    /// The exact two-sided sign test of `better` against `worse` (ties dropped).
    pub(crate) fn sign_p(&self) -> f64 {
        sign_test(self.better, self.worse)
    }
}

/// The exact two-sided binomial sign test: `P(X ≤ k) + P(X ≥ n − k)` for `X ~ B(n, ½)`,
/// `k = min(better, worse)`, `n = better + worse`, capped at 1; `1` for `n = 0`.
///
/// Up to `n = 1000` the terms are formed from `2⁻ⁿ`, which is exact, so a stratum's few records
/// give the closed form to the bit. Past that `2⁻ⁿ` nears binary64's underflow, and a pooled `n`
/// in the tens of thousands is summed in the log domain instead, good to about `n u` relative; a
/// `p` below binary64's range is returned as 0.
pub(crate) fn sign_test(better: usize, worse: usize) -> f64 {
    let n = better + worse;
    let k = better.min(worse);
    if n == 0 {
        return 1.0;
    }
    // Two-sided: the tails are equal by symmetry, and at `k = n/2` they overlap in one term,
    // which is why the cap is needed.
    let tail = match n <= DIRECT {
        true => tail_direct(n, k),
        false => tail_log(n, k),
    };
    (2.0 * tail).min(1.0)
}

/// The largest `n` whose `2⁻ⁿ` is a normal binary64 with room for the terms' growth.
const DIRECT: usize = 1000;

/// `P(X ≤ k)`, each term from the last by `C(n, i + 1) = C(n, i) (n − i) / (i + 1)`.
fn tail_direct(n: usize, k: usize) -> f64 {
    let mut term = 1.0f64;
    for _ in 0..n {
        term *= 0.5;
    }
    let mut tail = 0.0;
    for i in 0..=k {
        tail += term;
        term = term * (n - i) as f64 / (i + 1) as f64;
    }
    tail
}

/// `P(X ≤ k)` from `ln C(n, i) − n ln 2`, combined about the largest term.
fn tail_log(n: usize, k: usize) -> f64 {
    let ln2 = core::f64::consts::LN_2;
    let mut terms = Vec::with_capacity(k + 1);
    let mut ln_c = 0.0f64;
    for i in 0..=k {
        terms.push(ln_c - n as f64 * ln2);
        ln_c += ((n - i) as f64).ln() - ((i + 1) as f64).ln();
    }
    let top = terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    (top + terms.iter().map(|t| (t - top).exp()).sum::<f64>().ln()).exp()
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    fn tally(pairs: &[(Score, Score)]) -> Tally {
        let mut t = Tally::default();
        for &(a, b) in pairs {
            t.add(a, b);
        }
        t
    }

    fn f(x: f64) -> Score {
        Score::Finite(x)
    }

    #[test]
    fn the_sign_test_is_the_exact_binomial() {
        // Small n against the closed forms, 2 · C(n, ≤k) / 2ⁿ, to the bit.
        assert_eq!(sign_test(0, 0), 1.0);
        assert_eq!(sign_test(6, 0), 2.0 / 64.0);
        assert_eq!(sign_test(0, 6), 2.0 / 64.0);
        assert_eq!(sign_test(5, 1), 2.0 * 7.0 / 64.0);
        assert_eq!(sign_test(8, 2), 2.0 * 56.0 / 1024.0);
        // A balanced draw is 1, and a near-balanced one is capped there, not above.
        assert_eq!(sign_test(3, 3), 1.0);
        assert_eq!(sign_test(4, 3), 1.0);
        // `0038`'s published number: 104 worse of 125 moved rows is p = 2.1e-14.
        let p = sign_test(21, 104);
        assert!((2.05e-14..2.15e-14).contains(&p), "{p:e}");
        // Pooled sizes far past `C(n, k)`'s range stay finite, and a p past binary64's is 0.
        let p = sign_test(4_000, 6_000);
        assert!(p > 0.0 && p < 1e-80, "{p:e}");
        assert_eq!(sign_test(0, 5_000), 0.0);
        assert!((sign_test(50_000, 50_000) - 1.0).abs() < 1e-12);
        // The two summations agree where they meet.
        for k in [0, 10, 300, 450, 500] {
            let (d, l) = (tail_direct(DIRECT, k), tail_log(DIRECT, k));
            assert!((d - l).abs() <= 1e-11 * d, "{k}: {d:e} {l:e}");
        }
    }

    #[test]
    fn the_geometric_mean_is_of_the_records_that_differ() {
        let t = tally(&[(f(1.0), f(2.0)), (f(4.0), f(2.0)), (f(1.0), f(4.0))]);
        // (2 · ½ · 4)^(1/3) = 4^(1/3).
        let g = t.geomean().unwrap_or(f64::NAN);
        assert!((g - 4f64.cbrt()).abs() < 1e-15, "{g}");
        assert_eq!((t.better, t.worse, t.ties), (1, 2, 0));
        assert_eq!((t.worst_b, t.worst_a), (Some(4.0), Some(2.0)));
        // A tie moves no statistic but the count.
        let mut tied = t.clone();
        tied.add(f(3.0), f(3.0));
        assert_eq!(tied.geomean(), t.geomean());
        assert_eq!((tied.ties, tied.sign_p()), (1, t.sign_p()));
        assert_eq!(Tally::default().geomean(), None);
    }

    #[test]
    fn a_zero_error_is_a_direction_and_not_a_ratio() {
        let t = tally(&[(f(0.0), f(0.0)), (f(0.0), f(1.0)), (f(2.0), f(0.0))]);
        assert_eq!((t.ties, t.worse, t.better, t.zero_sided), (1, 1, 1, 2));
        assert_eq!(t.geomean(), None);
        assert_eq!(
            (t.worst_b, t.worst_a),
            (Some(f64::INFINITY), Some(f64::INFINITY))
        );
        assert_eq!((t.max_a, t.max_b), (Some(2.0), Some(1.0)));
        // Both maxima zero is a tie of the maxima, not 0/0.
        assert_eq!(tally(&[(f(0.0), f(0.0))]).max_ratio(), Some(1.0));
        assert_eq!(tally(&[(f(0.0), f(1.0))]).max_ratio(), Some(f64::INFINITY));
    }

    #[test]
    fn a_nonfinite_record_is_counted_and_compared_with_nothing() {
        let nan = Score::Finite(f64::NAN);
        let t = tally(&[
            (Score::NonFinite, f(1.0)),
            (f(1.0), Score::NonFinite),
            (Score::NonFinite, Score::NonFinite),
            (nan, f(1.0)),
        ]);
        assert_eq!((t.n, t.nonfinite_a, t.nonfinite_b), (4, 3, 2));
        assert_eq!((t.better, t.worse, t.ties), (0, 0, 0));
        assert_eq!((t.max_a, t.max_b, t.max_ratio()), (None, None, None));
        assert_eq!((t.geomean(), t.sign_p()), (None, 1.0));
    }

    #[test]
    fn pooling_adds_counts_and_keeps_the_extremes() {
        let mut s1 = tally(&[(f(1.0), f(2.0)), (f(1.0), f(2.0))]);
        let mut s2 = tally(&[(f(8.0), f(1.0)), (f(1.0), f(1.0))]);
        s1.close_stratum();
        s2.close_stratum();
        assert_eq!((s1.strata_worse, s2.strata_better), (1, 1));
        let mut id = Tally::default();
        id.merge(&s1);
        id.merge(&s2);
        assert_eq!((id.n, id.better, id.worse, id.ties), (4, 1, 2, 1));
        assert_eq!(
            (id.strata_better, id.strata_worse, id.strata_tied),
            (1, 1, 0)
        );
        assert_eq!(id.strata_ratio, Some((0.125, 2.0)));
        assert_eq!((id.worst_b, id.worst_a), (Some(2.0), Some(8.0)));
        // (2 · 2 · ⅛)^(1/3) = ½^(1/3).
        let g = id.geomean().unwrap_or(f64::NAN);
        assert!((g - 0.5f64.cbrt()).abs() < 1e-15, "{g}");
        assert_eq!(id.sign_p(), sign_test(1, 2));
    }
}
