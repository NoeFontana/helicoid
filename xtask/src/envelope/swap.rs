//! The two rules of the domination bar (`0088`): the paired swap test, whether a candidate's
//! maximum over a stratum exceeds the best oracle's by more than two equally good programs would,
//! and the gross-loss bound, which fails an excess no test on the maximum can resolve.
//!
//! The statistic is the bar's own, `T = max c − max o`, over the records both scored. Under the
//! null hypothesis that the two programs are equally good on the stratum, the two errors of each
//! record are exchangeable, so swapping them changes nothing; the p-value is `P(T' ≥ T)` over the
//! `2^n` swaps, each equally likely. Nothing is averaged. It is computed exactly, in closed form,
//! with no draw: a function of the errors alone, as everything else under D16 is.
//!
//! **The closed form.** Write each pair as its high `h_i = max(a_i, b_i)` and low `l_i`, sorted by
//! `h` descending. The side the first pair's high lands on has maximum `h_1`, each side with
//! probability ½. The other side's maximum `X` is `max(h_k, l_1, …, l_{k−1})` for the first `k ≥ 2`
//! whose high lands there, probability `2^-(k−1)`, and `max(l_1, …, l_n)` when none does; a pair
//! with `h = l` lands with probability 1, since both its errors are `h`. So
//! `T' = h_1 − X` or `X − h_1`, each with half of `X`'s probability, and `T' ≡ 0` when the first
//! pair is itself a tie. That is `O(n log n)`, against the `2^n` swaps it sums.

/// The level at which a loss is resolved (`0088` decision 1, the owner's answer to its question 2).
pub(super) const ALPHA: f64 = 0.001;

/// The gross-loss bound (`0088` decision 1): a maximum above `GROSS_FACTOR` times the best
/// oracle's plus [`GROSS_SLACK`] fails, resolved or not.
pub(super) const GROSS_FACTOR: f64 = 8.0;

/// The bound's absolute term, in `u`: an oracle that is exact on a stratum leaves room for one
/// ulp of rounding, and no more.
pub(super) const GROSS_SLACK: f64 = 1.0;

/// Whether a maximum `c` against the best oracle's `o` is a gross loss. The swap test cannot see a
/// loss confined to one record: every swap that moves that record lowers `T'`, every other keeps
/// it, so `p = ½` at any `n`. On today's corpus the largest ratio among the excesses it does not
/// resolve is 5.46 (`0088`), so 8 is outside the noise and a one-record catastrophe is not.
pub(super) fn gross(c: f64, o: f64) -> bool {
    c > GROSS_FACTOR * o + GROSS_SLACK
}

/// Whether `n` records can resolve anything when their errors were not read: no swap test on `n`
/// pairs gets below `2^-n`, so below 10 it never reaches [`ALPHA`].
pub(super) fn can_resolve(n: usize) -> bool {
    0.5_f64.powi(n as i32) <= ALPHA
}

/// The swap test's verdict on one stratum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Swap {
    /// `P(T' ≥ T)`.
    pub(super) p: f64,
    /// The least p-value any arrangement of these pairs could give: the probability of the largest
    /// `T'`. Above [`ALPHA`] the stratum cannot resolve a loss whatever its errors, so it is
    /// unresolvable; ties lift it, each tied pair being worth no swap.
    pub(super) floor: f64,
}

/// `T'` and its probability, over every swap of `pairs`, each `(candidate, oracle)` on one record.
fn distribution(pairs: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut sorted: Vec<(f64, f64)> = pairs.iter().map(|&(a, b)| (a.max(b), a.min(b))).collect();
    sorted.sort_by(|x, y| y.0.total_cmp(&x.0));
    let Some(&(h1, l1)) = sorted.first() else {
        return vec![(0.0, 1.0)];
    };
    // A tie within the first pair puts `h_1` on both sides: `T' = h_1 − h_1` at every swap.
    if h1.total_cmp(&l1).is_eq() {
        return vec![(h1 - l1, 1.0)];
    }
    let (mut other, mut low, mut mass) = (Vec::new(), l1, 1.0_f64);
    for &(h, l) in &sorted[1..] {
        if h.total_cmp(&l).is_eq() {
            other.push((h.max(low), mass));
            mass = 0.0;
            break;
        }
        mass *= 0.5;
        other.push((h.max(low), mass));
        low = low.max(l);
    }
    if mass > 0.0 {
        other.push((low, mass));
    }
    other
        .into_iter()
        .flat_map(|(x, q)| [(h1 - x, 0.5 * q), (x - h1, 0.5 * q)])
        .collect()
}

/// The swap test on `pairs`, each `(candidate, oracle)` on one record.
pub(super) fn test(pairs: &[(f64, f64)]) -> Swap {
    let max =
        |pick: fn(&(f64, f64)) -> f64| pairs.iter().map(pick).fold(f64::NEG_INFINITY, f64::max);
    let t = max(|p| p.0) - max(|p| p.1);
    let dist = distribution(pairs);
    let top = dist.iter().map(|d| d.0).fold(f64::NEG_INFINITY, f64::max);
    let at_least = |x: f64| dist.iter().filter(|d| d.0 >= x).map(|d| d.1).sum::<f64>();
    Swap {
        p: at_least(t),
        floor: at_least(top),
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    /// `n` values in `[1, 2)`, each the `splitmix64` hash of `(seed, index)`: independent
    /// streams for neighbouring seeds.
    fn errors(n: usize, seed: u64) -> Vec<f64> {
        let mix = |mut z: u64| {
            z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        };
        (0..n as u64)
            .map(|i| 1.0 + (mix(mix(seed) ^ i) >> 11) as f64 / (1_u64 << 53) as f64)
            .collect()
    }

    /// Every swap, enumerated: the definition the closed form is checked against.
    fn enumerated(pairs: &[(f64, f64)]) -> Swap {
        let stat = |m: u64| {
            let (mut c, mut o) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
            for (i, &(a, b)) in pairs.iter().enumerate() {
                let (a, b) = if m >> i & 1 == 1 { (b, a) } else { (a, b) };
                c = c.max(a);
                o = o.max(b);
            }
            c - o
        };
        let all = 1_u64 << pairs.len();
        let ts: Vec<f64> = (0..all).map(stat).collect();
        let top = ts.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let share = |x: f64| ts.iter().filter(|&&t| t >= x).count() as f64 / all as f64;
        Swap {
            p: share(stat(0)),
            floor: share(top),
        }
    }

    /// The closed form is the enumeration to the bit, for every `n` up to 12, on distinct errors,
    /// on errors with ties across pairs (small integers), and on pairs tied within themselves.
    #[test]
    fn the_closed_form_is_the_enumeration_to_the_bit() {
        for n in 1..=12 {
            for seed in 0..40_u64 {
                let (a, b) = (errors(n, 2 * seed + 1), errors(n, 2 * seed + 2));
                let distinct: Vec<(f64, f64)> = a.iter().copied().zip(b.iter().copied()).collect();
                let small: Vec<(f64, f64)> = distinct
                    .iter()
                    .map(|&(x, y)| ((x * 4.0).floor(), (y * 4.0).floor()))
                    .collect();
                let mut tied = distinct.clone();
                for (i, p) in tied.iter_mut().enumerate() {
                    if (i as u64 + seed) % 3 == 0 {
                        p.1 = p.0;
                    }
                }
                for pairs in [&distinct, &small, &tied] {
                    let (got, want) = (test(pairs), enumerated(pairs));
                    assert_eq!(got.p.to_bits(), want.p.to_bits(), "p, n = {n}: {pairs:?}");
                    assert_eq!(got.floor.to_bits(), want.floor.to_bits(), "floor, n = {n}");
                }
            }
        }
    }

    #[test]
    fn a_stratum_needs_ten_untied_records_to_resolve() {
        assert!(!can_resolve(4) && !can_resolve(9) && can_resolve(10) && can_resolve(801));
        let pairs: Vec<(f64, f64)> = errors(12, 3).into_iter().map(|e| (2.0 * e, e)).collect();
        assert_eq!(test(&pairs).floor, 2.0_f64.powi(-12));
        // Three pairs tied within themselves are worth no swap: the floor rises to `2^-9`, past
        // `ALPHA`.
        let mut tied = pairs;
        for p in &mut tied[..3] {
            p.1 = p.0 / 4.0;
            p.0 = p.1;
        }
        assert_eq!(test(&tied).floor, 2.0_f64.powi(-9));
        assert!(test(&tied).floor > ALPHA);
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
            let max = |f: fn(&(f64, f64)) -> f64| pairs.iter().map(f).fold(0.0, f64::max);
            if max(|p| p.0) > max(|p| p.1) {
                above += 1;
                failed += usize::from(test(&pairs).p <= ALPHA);
            }
        }
        assert!((60..140).contains(&above), "{above} of 200 above");
        assert_eq!(failed, 0);
    }

    /// A program worse on every record by twice: only the identity reaches `T`, so the p-value is
    /// the floor, `2^-n` exactly, resolved at 12 and 64 records and not at 6.
    #[test]
    fn a_loss_on_every_record_resolves_from_ten_records() {
        let pairs = |n: usize| -> Vec<(f64, f64)> {
            errors(n, 7).into_iter().map(|e| (2.0 * e, e)).collect()
        };
        assert_eq!(test(&pairs(64)).p, 2.0_f64.powi(-64));
        assert_eq!(test(&pairs(12)).p, 2.0_f64.powi(-12));
        assert_eq!(test(&pairs(6)).p, 2.0_f64.powi(-6));
    }

    /// A loss on one record is invisible to the test, `p = ½`, at any `n`; the gross-loss bound
    /// is what fails it, and it passes an equal program and every excess `0088` measured.
    #[test]
    fn a_one_record_catastrophe_is_the_bounds_to_fail() {
        for n in [12, 64, 801] {
            let mut pairs: Vec<(f64, f64)> = errors(n, 9).into_iter().map(|e| (e, e)).collect();
            pairs[n / 2].0 = 1e6;
            assert_eq!(test(&pairs).p, 0.5, "n = {n}");
            assert!(gross(1e6, 2.0));
        }
        assert!(!gross(2.0, 2.0) && !gross(6.66, 1.221) && !gross(0.9, 0.0));
        assert!(gross(1.5, 0.0) && gross(256.1, 2.132));
    }
}
