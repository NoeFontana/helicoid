//! The paired swap test of the domination bar (`0088`): whether a candidate's maximum over a
//! stratum exceeds the best oracle's by more than two equally good programs would.
//!
//! The statistic is the bar's own, `T = max c − max o`, over the records both scored. Under the
//! null hypothesis that the two programs are equally good on the stratum, the two errors of each
//! record are exchangeable, so swapping them changes nothing; the p-value is `P(T' ≥ T)` over the
//! swaps. Nothing is averaged. Up to [`EXACT_UP_TO`] records every swap is enumerated; above, [`DRAWS`]
//! swaps are drawn from `splitmix64` at [`SEED`], so the verdict is a function of the errors
//! alone, as everything else under D16 is.

/// The level at which a loss is resolved (`0088` decision 1, the owner's answer to its question 2).
pub(super) const ALPHA: f64 = 0.001;

/// Every swap is enumerated up to this many records: `2^16` of them.
pub(super) const EXACT_UP_TO: usize = 16;

/// The swaps drawn above [`EXACT_UP_TO`] records.
pub(super) const DRAWS: u64 = 1 << 16;

/// The generator's seed (`conformance/generate/gen/rng.py`, "helicoid").
pub(super) const SEED: u64 = 0x6865_6C69_636F_6964;

/// Whether a stratum of `n` records can resolve a loss at all: no swap test on `n` pairs reaches a
/// p-value below `2^-n`, so below 10 records it never reaches [`ALPHA`] (`0088` decision 2).
pub(super) fn resolvable(n: usize) -> bool {
    n >= 64 || 0.5_f64.powi(n as i32) <= ALPHA
}

/// `splitmix64`, the generator's stream.
struct SplitMix64(u64);

impl SplitMix64 {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

/// `max c' − max o'` with the records whose bit is set in `swapped` exchanged.
fn statistic(pairs: &[(f64, f64)], swapped: &dyn Fn(usize) -> bool) -> f64 {
    let (mut c, mut o) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    for (i, &(a, b)) in pairs.iter().enumerate() {
        let (a, b) = if swapped(i) { (b, a) } else { (a, b) };
        c = c.max(a);
        o = o.max(b);
    }
    c - o
}

/// The p-value of `pairs`, each `(candidate, oracle)` on one record: `P(T' ≥ T)` over the swaps,
/// the identity among them, so never below `2^-n` exactly or `1/(DRAWS + 1)` by draws.
pub(super) fn p_value(pairs: &[(f64, f64)]) -> f64 {
    let t = statistic(pairs, &|_| false);
    let n = pairs.len();
    if n <= EXACT_UP_TO {
        let all = 1_u64 << n;
        let at_least = (0..all)
            .filter(|&m| statistic(pairs, &|i| m >> i & 1 == 1) >= t)
            .count();
        return at_least as f64 / all as f64;
    }
    let mut rng = SplitMix64(SEED);
    let words = n.div_ceil(64);
    let mut bits = vec![0_u64; words];
    let mut at_least = 1_u64;
    for _ in 0..DRAWS {
        for w in &mut bits {
            *w = rng.next();
        }
        if statistic(pairs, &|i| bits[i / 64] >> (i % 64) & 1 == 1) >= t {
            at_least += 1;
        }
    }
    at_least as f64 / (DRAWS + 1) as f64
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    /// `n` errors in `[1, 2)` from a stream of their own: what two equally good programs read.
    fn errors(n: usize, seed: u64) -> Vec<f64> {
        let mut rng = SplitMix64(seed);
        (0..n)
            .map(|_| 1.0 + (rng.next() >> 11) as f64 / (1_u64 << 53) as f64)
            .collect()
    }

    #[test]
    fn the_stream_is_the_generators() {
        // `gen.rng.SplitMix64(SEED).next_u64()`, computed by the generator's own code.
        assert_eq!(SplitMix64(SEED).next(), 0x11b3_0f44_c266_74b9);
    }

    #[test]
    fn only_a_stratum_of_ten_records_or_more_can_resolve() {
        assert!(!resolvable(4) && !resolvable(6) && !resolvable(9));
        assert!(resolvable(10) && resolvable(12) && resolvable(64) && resolvable(801));
    }

    /// Exchangeable errors never resolve at `ALPHA`: over 200 pairs of equally good programs, the
    /// candidate's maximum is above half the time, and the test fails none of them.
    #[test]
    fn two_equally_good_programs_never_fail() {
        let (mut above, mut failed) = (0, 0);
        for k in 0..200_u64 {
            let n = [12, 24, 64][(k % 3) as usize];
            let (c, o) = (errors(n, 2 * k + 1), errors(n, 2 * k + 2));
            let pairs: Vec<(f64, f64)> = c.into_iter().zip(o).collect();
            if statistic(&pairs, &|_| false) > 0.0 {
                above += 1;
                failed += usize::from(p_value(&pairs) <= ALPHA);
            }
        }
        assert!((60..140).contains(&above), "{above} of 200 above");
        assert_eq!(failed, 0);
    }

    /// A program worse on every record by twice: only the identity reaches `T`, since any swap
    /// puts an error of at least 2 among the oracle's, all below 2. So the p-value is the floor,
    /// `2^-n` exactly, resolved at 12 and 64 records and not at 6.
    #[test]
    fn a_loss_on_every_record_resolves_from_ten_records() {
        let pairs = |n: usize| -> Vec<(f64, f64)> {
            errors(n, 7).into_iter().map(|e| (2.0 * e, e)).collect()
        };
        assert_eq!(p_value(&pairs(64)), 1.0 / (DRAWS + 1) as f64);
        assert_eq!(p_value(&pairs(12)), 2.0_f64.powi(-12));
        assert_eq!(p_value(&pairs(6)), 2.0_f64.powi(-6));
        assert!(p_value(&pairs(12)) <= ALPHA && p_value(&pairs(6)) > ALPHA);
    }

    /// The draws are a function of the errors: two runs give the same p-value to the bit.
    #[test]
    fn the_drawn_p_value_is_deterministic() {
        let c = errors(40, 11);
        let o = errors(40, 12);
        let pairs: Vec<(f64, f64)> = c.into_iter().zip(o).collect();
        assert_eq!(p_value(&pairs).to_bits(), p_value(&pairs).to_bits());
    }
}
