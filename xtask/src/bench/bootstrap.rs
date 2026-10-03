//! The paired bootstrap 95% CI of a ratio (`docs/PHASE1.md` §9).
//!
//! Criterion takes 100 samples on a linear schedule, so sample `i` of a run is the `i`-th block of
//! its measurement window. Two runs pair **by position in that schedule**, not by iteration count:
//! the counts differ run to run, because criterion sizes them from what warm-up measured. Position
//! is the right pairing even so — it matches block `i` with block `i`, which is what carries the
//! within-run trend (frequency ramp, cache state), and an independent resample would charge that
//! trend to the interval as if it were noise.
//!
//! **The statistic is the ratio of the runs' median per-iteration times, not their means.** §9
//! prescribes "the paired bootstrap 95% CI of the ratio" and does not name the location statistic,
//! so it is chosen by measurement. A window's samples are the true cost plus one-sided
//! contamination — a preempted sample costs more, never less — so the mean is not an estimator of
//! the cost. Measured on this host, the same benchmark run eight times back to back:
//!
//! | statistic of a window | spread across the eight windows |
//! |---|---|
//! | mean | 2.97% |
//! | **median** | **1.04%** |
//! | 10th percentile | 0.91% |
//! | minimum | 1.41% |
//!
//! Within a single window the interquartile spread was 1.1% to 2.7% while the extremes spread up
//! to 96%: a handful of samples at 1.3x to 2x the typical cost, and the rest extremely tight. The
//! mean carries those; worse, resampling them with replacement lets their multiplicity swing every
//! resample, which is how a benchmark whose true cost is reproducible to 1% earned an A/A floor of
//! 38%. The median takes the tight middle and has a 50% breakdown point. The 10th percentile scored
//! marginally better and the minimum worse — the minimum is one sample, so it is the least stable
//! of the robust three.
//!
//! What is required of the two runs is the same *number* of samples; an unequal count means the two
//! did not measure the same way and there is no pairing to resample.
//!
//! Resampling is seeded (splitmix64, the corpus generator's), so a verdict is a function of the
//! samples and the seed alone and two runs of the gate on one pair of files agree to the bit.

/// splitmix64, as `conformance/generate` uses: a seeded stream with no dependency.
pub(super) struct Rng(u64);

impl Rng {
    pub(super) fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform on `0..n`, by Lemire's multiply-shift; `n > 0`.
    fn below(&mut self, n: usize) -> usize {
        debug_assert!(n > 0, "Rng::below: empty range");
        // The bias of a single multiply is at most `n / 2^64`, far below the 2^-20 resolution of
        // a 1024-resample percentile, so no rejection loop.
        ((u128::from(self.next_u64()) * n as u128) >> 64) as usize
    }
}

/// The median of a whole sample, by the same nearest rank `ratio_ci` takes.
///
/// A window's cost, for choosing between the alternations of one benchmark: the least costly pair
/// is the least disturbed, because contamination is one-sided.
pub(super) fn median(per_iter: &[f64]) -> f64 {
    let idx: Vec<usize> = (0..per_iter.len()).collect();
    let mut buf = Vec::with_capacity(per_iter.len());
    median_at(per_iter, &idx, &mut buf)
}

/// A 95% confidence interval of a ratio.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Interval {
    pub(super) low: f64,
    pub(super) point: f64,
    pub(super) high: f64,
}

/// How many resamples a CI is taken over. 1024 puts the 2.5% and 97.5% ranks at 25 and 998, both
/// interior, so neither bound is an extreme order statistic.
pub(super) const RESAMPLES: usize = 1024;

/// The median of `times[i] / iters[i]` over the given indices, by nearest rank.
///
/// `buf` is the caller's scratch so that 1024 resamples do not allocate 1024 times. Nearest rank
/// and no interpolation, the convention `ratio_ci`'s percentiles already take.
fn median_at(per_iter: &[f64], idx: &[usize], buf: &mut Vec<f64>) -> f64 {
    buf.clear();
    buf.extend(idx.iter().map(|&i| per_iter[i]));
    buf.sort_by(f64::total_cmp);
    let rank = ((0.5 * buf.len() as f64).ceil() as usize).clamp(1, buf.len()) - 1;
    buf.get(rank).copied().unwrap_or(f64::NAN)
}

/// The paired bootstrap CI of `candidate / baseline`, both as per-iteration times in sample order.
///
/// `point` is the ratio of the full-sample medians. Returns `None` when the two differ in length
/// (the schedules did not match, so the samples are not pairs) or either is empty.
pub(super) fn ratio_ci(baseline: &[f64], candidate: &[f64], seed: u64) -> Option<Interval> {
    if baseline.len() != candidate.len() || baseline.is_empty() {
        return None;
    }
    let n = baseline.len();
    let all: Vec<usize> = (0..n).collect();
    let mut buf = Vec::with_capacity(n);
    let point = median_at(candidate, &all, &mut buf) / median_at(baseline, &all, &mut buf);

    let mut rng = Rng::new(seed);
    let mut ratios = Vec::with_capacity(RESAMPLES);
    let mut idx = vec![0usize; n];
    for _ in 0..RESAMPLES {
        for slot in idx.iter_mut() {
            *slot = rng.below(n);
        }
        let b = median_at(baseline, &idx, &mut buf);
        if b > 0.0 {
            ratios.push(median_at(candidate, &idx, &mut buf) / b);
        }
    }
    if ratios.len() < RESAMPLES {
        return None;
    }
    ratios.sort_by(f64::total_cmp);
    // Nearest-rank percentiles, as `error-analysis.md` EA.13(d) takes them: no interpolation.
    let rank = |p: f64| ratios[((p * RESAMPLES as f64).ceil() as usize).clamp(1, RESAMPLES) - 1];
    Some(Interval {
        low: rank(0.025),
        point,
        high: rank(0.975),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The stream is the generator's: a fixed seed gives a fixed sequence, so a verdict is
    /// reproducible.
    #[test]
    fn the_rng_is_seeded_and_bounded() {
        let mut a = Rng::new(7);
        let mut b = Rng::new(7);
        assert_eq!(a.next_u64(), b.next_u64());
        let mut r = Rng::new(1);
        assert!((0..1000).all(|_| r.below(5) < 5));
        assert_eq!(Rng::new(1).below(1), 0);
    }

    fn ci(b: &[f64], c: &[f64], seed: u64) -> Result<Interval, String> {
        ratio_ci(b, c, seed).ok_or_else(|| "not a pair".to_string())
    }

    /// Identical inputs: the ratio is exactly 1 and every resample is 1, so the CI is degenerate.
    #[test]
    fn identical_samples_give_a_unit_interval() -> Result<(), String> {
        let s: Vec<f64> = (1..=100).map(f64::from).collect();
        let i = ci(&s, &s, 1)?;
        for x in [i.low, i.point, i.high] {
            assert_eq!(x.to_bits(), 1.0_f64.to_bits(), "{i:?}");
        }
        Ok(())
    }

    /// A uniform scaling is recovered exactly, bounds included: every resample sees the same ratio.
    #[test]
    fn a_uniform_scaling_is_the_whole_interval() -> Result<(), String> {
        let b: Vec<f64> = (1..=100).map(f64::from).collect();
        let c: Vec<f64> = b.iter().map(|x| x * 1.25).collect();
        let i = ci(&b, &c, 1)?;
        for x in [i.low, i.point, i.high] {
            assert!((x - 1.25).abs() < 1e-12, "{i:?}");
        }
        Ok(())
    }

    /// Noise widens the interval around 1 and the point estimate stays near it.
    #[test]
    fn noise_widens_the_interval_around_one() -> Result<(), String> {
        let mut r = Rng::new(42);
        let jitter = |r: &mut Rng| 1.0 + (r.below(2001) as f64 - 1000.0) / 10_000.0;
        let b: Vec<f64> = (0..100).map(|_| 10.0 * jitter(&mut r)).collect();
        let c: Vec<f64> = (0..100).map(|_| 10.0 * jitter(&mut r)).collect();
        let i = ci(&b, &c, 1)?;
        assert!(i.low < 1.0 && i.high > 1.0, "{i:?}");
        assert!((i.point - 1.0).abs() < 0.05, "{i:?}");
        Ok(())
    }

    /// The load-bearing property of the median statistic: one-sided contamination of the
    /// candidate's window — a few preempted samples, which is what the host actually produces —
    /// does not move the verdict, where the mean would report a slowdown that is not there.
    #[test]
    fn a_contaminated_window_does_not_move_the_median_ratio() -> Result<(), String> {
        let b: Vec<f64> = (0..100).map(|i| 10.0 + f64::from(i % 3) / 100.0).collect();
        let mut c = b.clone();
        // Eight samples at double cost, the shape measured on this host (extremes 1.3x to 2x).
        for x in c.iter_mut().take(8) {
            *x *= 2.0;
        }
        let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
        // The mean of the same samples reports an 8% slowdown.
        assert!((mean(&c) / mean(&b) - 1.08).abs() < 0.01);
        let i = ci(&b, &c, 1)?;
        assert!(i.high < 1.01, "{i:?}");
        Ok(())
    }

    /// Unequal lengths are not pairs: the schedules differed, so there is no verdict to give.
    #[test]
    fn unpaired_or_empty_samples_have_no_interval() {
        assert!(ratio_ci(&[1.0, 2.0], &[1.0], 1).is_none());
        assert!(ratio_ci(&[], &[], 1).is_none());
    }

    /// A seed changes the bounds and not the point estimate, which is the full sample's.
    #[test]
    fn the_point_estimate_does_not_depend_on_the_seed() -> Result<(), String> {
        let mut r = Rng::new(9);
        let b: Vec<f64> = (0..100).map(|_| 5.0 + r.below(100) as f64 / 50.0).collect();
        let c: Vec<f64> = (0..100).map(|_| 5.0 + r.below(100) as f64 / 50.0).collect();
        let (x, y) = (ci(&b, &c, 1)?, ci(&b, &c, 2)?);
        assert_eq!(x.point.to_bits(), y.point.to_bits());
        Ok(())
    }
}
