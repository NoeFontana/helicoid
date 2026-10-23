//! `helicoid:portable-atan2`: the shipped program at [`Portable`], whose `atan2` is the portable
//! candidate of `0105` (draft) and whose every other operation is the shipped scalar's. The kernels
//! and their fit are [`kernel`]'s.
//!
//! Planted, as `helicoid:host-std` is: a twin is not a candidate, and no bar reads its rows. It
//! answers both precisions, since every operation but `atan2` is the shipped scalar's own.

pub(crate) mod kernel;

use helicoid_linalg::Precision;

use super::{program, Helicoid};
use crate::conformance::corpus::Record;
use crate::conformance::subject::{Output, Registered, Subject};
use kernel::{fits, Portable};

/// The subject's name, and the value of `--subject`.
pub(crate) const NAME: &str = "helicoid:portable-atan2";

/// The shipped program at [`Portable`].
pub(crate) struct Twin;

impl Subject for Twin {
    fn name(&self) -> &str {
        NAME
    }

    fn supports(&self, fn_id: &str) -> bool {
        Helicoid.supports(fn_id)
    }

    fn eval(&self, fn_id: &str, record: &Record, precision: Precision) -> Output {
        if fits().is_err() {
            return Output::new();
        }
        program::<Portable<f64>, Portable<f32>>(fn_id, record, precision)
    }
}

pub(crate) fn registered() -> Registered {
    let version = format!("{}@portable-atan2", env!("CARGO_PKG_VERSION"));
    Registered {
        version: version.clone(),
        version_f32: version,
        no_f32: None,
        planted: true,
        subject: Box::new(Twin),
    }
}

#[cfg(test)]
mod tests {
    use super::kernel::{atan2_f32, atan2_f64, centre, fits, CENTRES};

    /// splitmix64, so a draw is the same on every host.
    struct Draw(u64);

    impl Draw {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }

        /// `±m·2^e`, `m` uniform in `[1, 2)`, `e` uniform in `[lo, hi]`.
        fn wide(&mut self, lo: i32, hi: i32) -> f64 {
            let r = self.next();
            let span = u64::try_from(hi - lo + 1).unwrap_or(1);
            let e = lo + i32::try_from(r % span).unwrap_or(0);
            let m = 1.0 + (r >> 11) as f64 / (1u64 << 53) as f64;
            let x = libm::ldexp(m, e);
            if r & 1 == 0 {
                x
            } else {
                -x
            }
        }
    }

    /// A double-double reference, about `2^-100` relative: enough to read an error in binary64
    /// ulps to many digits. `π/4` is its own `atan(1)`, so nothing here is a typed constant.
    mod dd {
        #[derive(Clone, Copy, Debug)]
        pub(super) struct Dd(pub(super) f64, pub(super) f64);

        const SPLIT: f64 = 134_217_729.0;

        fn fast_two_sum(a: f64, b: f64) -> (f64, f64) {
            let s = a + b;
            (s, b - (s - a))
        }

        /// Dekker's exact product, for operands whose halves' products are normal.
        fn two_product(a: f64, b: f64) -> (f64, f64) {
            let halves = |x: f64| {
                let c = SPLIT * x;
                let hi = c - (c - x);
                (hi, x - hi)
            };
            let h = a * b;
            let ((ah, al), (bh, bl)) = (halves(a), halves(b));
            (h, ((ah * bh - h) + ah * bl + al * bh) + al * bl)
        }

        fn two_sum(a: f64, b: f64) -> (f64, f64) {
            let s = a + b;
            let bb = s - a;
            (s, (a - (s - bb)) + (b - bb))
        }

        impl Dd {
            pub(super) fn of(x: f64) -> Self {
                Dd(x, 0.0)
            }
            pub(super) fn add(self, o: Self) -> Self {
                let (s, e) = two_sum(self.0, o.0);
                let (t, f) = two_sum(self.1, o.1);
                let (s, e) = fast_two_sum(s, e + t);
                let (s, e) = fast_two_sum(s, e + f);
                Dd(s, e)
            }
            pub(super) fn neg(self) -> Self {
                Dd(-self.0, -self.1)
            }
            pub(super) fn mul(self, o: Self) -> Self {
                let (p, e) = two_product(self.0, o.0);
                let (p, e) = fast_two_sum(p, e + (self.0 * o.1 + self.1 * o.0));
                Dd(p, e)
            }
            pub(super) fn div(self, o: Self) -> Self {
                let q1 = self.0 / o.0;
                let r = self.add(o.mul(Dd::of(q1)).neg());
                let q2 = r.0 / o.0;
                let r = r.add(o.mul(Dd::of(q2)).neg());
                let q3 = r.0 / o.0;
                Dd::of(q1).add(Dd::of(q2)).add(Dd::of(q3))
            }
            pub(super) fn sqrt(self) -> Self {
                let s = self.0.sqrt();
                let r = self.add(Dd::of(s).mul(Dd::of(s)).neg());
                Dd::of(s).add(Dd::of(r.0 / (2.0 * s)))
            }
        }

        /// `atan(t)` for `0 ≤ t ≤ 1`: four half-angle steps, then the series.
        pub(super) fn atan(t: Dd) -> Dd {
            let one = Dd::of(1.0);
            let mut z = t;
            for _ in 0..4 {
                z = z.div(one.add(one.add(z.mul(z)).sqrt()));
            }
            let z2 = z.mul(z);
            let (mut sum, mut pow) = (Dd::of(0.0), z);
            for j in 0..20 {
                let term = pow.div(Dd::of(f64::from(2 * j + 1)));
                sum = if j % 2 == 0 {
                    sum.add(term)
                } else {
                    sum.add(term.neg())
                };
                pow = pow.mul(z2);
            }
            sum.mul(Dd::of(16.0))
        }

        /// `atan2(y, x)` for finite, non-zero arguments.
        pub(super) fn atan2(y: f64, x: f64) -> Dd {
            let (ax, ay) = (x.abs(), y.abs());
            let (n, d) = if ay > ax { (ax, ay) } else { (ay, ax) };
            let a = atan(Dd::of(n).div(Dd::of(d)));
            let pi4 = atan(Dd::of(1.0));
            let pi2 = pi4.mul(Dd::of(2.0));
            let pi = pi4.mul(Dd::of(4.0));
            let a = match (x < 0.0, ay > ax) {
                (false, false) => a,
                (false, true) => pi2.add(a.neg()),
                (true, true) => pi2.add(a),
                (true, false) => pi.add(a.neg()),
            };
            if y < 0.0 {
                a.neg()
            } else {
                a
            }
        }
    }

    /// `|got − want|` in ulps of `want` at `bits` significand bits.
    fn ulps(got: f64, want: dd::Dd, bits: i32) -> f64 {
        let e = libm::ilogb(want.0).max(if bits == 53 { -1022 } else { -126 });
        ((got - want.0) - want.1).abs() / libm::ldexp(1.0, e - (bits - 1))
    }

    /// The worst error of `subject` at binary64 over `draws` pairs, and where.
    fn f64_error(
        draws: usize,
        seed: u64,
        (lo, hi): (i32, i32),
        subject: impl Fn(f64, f64) -> Result<f64, String>,
    ) -> Result<(f64, f64, f64), String> {
        let mut draw = Draw(seed);
        let mut worst = (0.0, 0.0, 0.0);
        for _ in 0..draws {
            let (y, x) = (draw.wide(lo, hi), draw.wide(lo, hi));
            let err = ulps(subject(y, x)?, dd::atan2(y, x), 53);
            if err > worst.0 {
                worst = (err, y, x);
            }
        }
        Ok(worst)
    }

    /// The tuning readout: the candidate's, `libm`'s and the host's worst binary64 error against
    /// the double-double reference, near and far ratios.
    ///
    /// `cargo nextest run -p xtask --run-ignored only --no-capture -- portable_atan2`.
    #[test]
    #[ignore = "a measurement, printed for 0105 (draft) to cite"]
    #[allow(clippy::print_stdout)]
    fn the_binary64_kernel_reads_its_error() -> Result<(), String> {
        for (range, label) in [((-2, 2), "|y/x| in 2^±4"), ((-30, 30), "|y/x| in 2^±60")] {
            let (err, y, x) = f64_error(1 << 20, 5, range, at64)?;
            println!("f64 {label}: max {err:.4} ulp at atan2({y:e}, {x:e})");
            let (err, y, x) = f64_error(1 << 20, 5, range, |y, x| Ok(libm::atan2(y, x)))?;
            println!("libm::atan2 {label}: max {err:.4} ulp at atan2({y:e}, {x:e})");
            let (err, y, x) = f64_error(1 << 20, 5, range, |y, x| Ok(y.atan2(x)))?;
            println!("host std {label}: max {err:.4} ulp at atan2({y:e}, {x:e})");
        }
        Ok(())
    }

    /// The candidate's worst binary64 error per centre `j` and quadrant `q`.
    #[test]
    #[ignore = "a measurement, printed for 0105 (draft) to cite"]
    #[allow(clippy::print_stdout)]
    fn the_binary64_kernel_reads_its_error_per_arm() -> Result<(), String> {
        let mut draw = Draw(6);
        let mut worst = std::collections::BTreeMap::<(u8, u8), (f64, f64)>::new();
        for _ in 0..1 << 20 {
            let (y, x) = (draw.wide(-2, 2), draw.wide(-2, 2));
            let (n, d) = (y.abs().min(x.abs()), y.abs().max(x.abs()));
            let (c, ..) = centre(n, d);
            let j = (0u8..=8)
                .find(|&j| (f64::from(j) / 8.0).to_bits() == c.to_bits())
                .unwrap_or(9);
            let q = match (y.abs() > x.abs(), x < 0.0) {
                (true, _) => 2,
                (false, true) => 4,
                (false, false) => 0,
            };
            let err = ulps(at64(y, x)?, dd::atan2(y, x), 53);
            let cell = worst.entry((j, q)).or_insert((0.0, 0.0));
            if err > cell.0 {
                *cell = (err, n / d);
            }
        }
        for ((j, q), (err, t)) in &worst {
            println!("j {j} q {q}: max {err:.4} ulp at n/d = {t:.4}");
        }
        Ok(())
    }

    /// The reference against itself where the answer is known: `atan2(1, 1) = π/4` and
    /// `atan2(1, 0+) ≈ π/2` both round to `core`'s constants.
    #[test]
    fn the_reference_rounds_to_the_known_constants() {
        let pi4 = dd::atan2(1.0, 1.0);
        assert_eq!(pi4.0.to_bits(), core::f64::consts::FRAC_PI_4.to_bits());
        assert!(ulps(core::f64::consts::FRAC_PI_4, pi4, 53) < 0.5);
        let pi = dd::atan2(1e-300, -1.0);
        assert_eq!(pi.0.to_bits(), core::f64::consts::PI.to_bits());
    }

    fn at64(y: f64, x: f64) -> Result<f64, String> {
        Ok(atan2_f64(y, x, &fits()?.f64))
    }

    fn at32(y: f32, x: f32) -> Result<f32, String> {
        let f = fits()?;
        Ok(atan2_f32(y, x, &f.f32, (f.f64.pi4_hi, f.f64.pi4_lo)))
    }

    fn same64(a: f64, b: f64) -> bool {
        a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan())
    }

    fn same32(a: f32, b: f32) -> bool {
        a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan())
    }

    /// The last centre is `1`, so its `atan` is the split of `π/4`; `q·π/4_hi` is exact for
    /// `q ≤ 4` because `π/4_hi` has two spare low bits.
    #[test]
    fn the_committed_fit_ends_at_pi_over_4_with_two_spare_bits() -> Result<(), String> {
        let f = &fits()?.f64;
        assert_eq!(f.atan_hi[CENTRES].to_bits(), f.pi4_hi.to_bits());
        assert_eq!(f.atan_lo[CENTRES].to_bits(), f.pi4_lo.to_bits());
        assert!(f.pi4_hi.to_bits().trailing_zeros() >= 2);
        Ok(())
    }

    /// `centre`'s claims, against the double-double reference: `num = n − c·d` exactly, and
    /// `den + den_lo = d + c·n` to `2⁻¹⁰⁴`, on draws across the window and at every midpoint.
    #[test]
    fn the_centred_numerator_is_exact() {
        let mut draw = Draw(8);
        for i in 0u32..200_000 {
            let d = draw.wide(-460, 989).abs();
            let n = if i % 2 == 0 {
                d * draw.wide(-60, -1).abs()
            } else {
                let mid = f64::from(2 * (i / 2 % 8) + 1) / 16.0;
                let k = draw.next() % 5;
                f64::from_bits(((d * mid).to_bits() + k) - 2)
            }
            .min(d);
            let (c, num, den, den_lo) = centre(n, d);
            assert!(
                (n / d - c).abs() <= 1.0 / 16.0 + 1e-15,
                "n {n:e} d {d:e} c {c}"
            );
            let exact = dd::Dd::of(n).add(dd::Dd::of(c).mul(dd::Dd::of(d)).neg());
            assert_eq!(
                (num.to_bits(), exact.1),
                (exact.0.to_bits(), 0.0),
                "n {n:e} d {d:e}"
            );
            let want = dd::Dd::of(d).add(dd::Dd::of(c).mul(dd::Dd::of(n)));
            let got = dd::Dd::of(den).add(dd::Dd::of(den_lo));
            let off = got.add(want.neg());
            assert!(
                (off.0 + off.1).abs() <= libm::ldexp(den, -104),
                "n {n:e} d {d:e}"
            );
        }
    }

    /// Wherever an argument is zero, infinite or NaN, the answer is IEEE 754's, correctly rounded:
    /// `±0`, `±π/2`, `±π`, `±π/4`, `±3π/4`, or NaN. At `f64` that is `libm`'s to the bit. At `f32`
    /// it is `libm`'s binary64 answer rounded, since `libm::atan2f` returns `π/2` one ulp low where
    /// `y` is infinite or dwarfs `x` (`0x3fc90fda` against the nearest, `0x3fc90fdb`); those
    /// binary64 constants are nowhere near a binary32 midpoint.
    #[test]
    fn special_values_are_ieee_s() -> Result<(), String> {
        let mut atan2f_off = Vec::new();
        let edge = [0.0, -0.0, f64::INFINITY, f64::NEG_INFINITY, f64::NAN];
        let any = [
            0.0,
            -0.0,
            1.0,
            -1.0,
            3.5e-310,
            -2.0e300,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NAN,
        ];
        for &a in &edge {
            for &b in &any {
                for (y, x) in [(a, b), (b, a)] {
                    let (got, want) = (at64(y, x)?, libm::atan2(y, x));
                    assert!(
                        same64(got, want),
                        "f64 atan2({y:e}, {x:e}) = {got:e}, libm {want:e}"
                    );
                    let (y, x) = (y as f32, x as f32);
                    let want = libm::atan2(f64::from(y), f64::from(x)) as f32;
                    let got = at32(y, x)?;
                    assert!(
                        same32(got, want),
                        "f32 atan2({y:e}, {x:e}) = {got:e}, want {want:e}"
                    );
                    if !same32(libm::atan2f(y, x), want) {
                        atan2f_off.push((y, x));
                    }
                }
            }
        }
        let inf = f32::INFINITY;
        assert_eq!(atan2f_off, [(inf, 1.0), (-inf, 1.0)]);
        Ok(())
    }

    /// `atan2(−y, x) = −atan2(y, x)` bit for bit: the kernel works on `|y|` and restores its sign.
    #[test]
    fn atan2_is_odd_in_y() -> Result<(), String> {
        let mut draw = Draw(1);
        for _ in 0..100_000 {
            let (y, x) = (draw.wide(-1074, 1023), draw.wide(-1074, 1023));
            assert!(same64(at64(-y, x)?, -at64(y, x)?), "({y:e}, {x:e})");
            let (y, x) = (y as f32, x as f32);
            assert!(same32(at32(-y, x)?, -at32(y, x)?), "({y:e}, {x:e})");
        }
        Ok(())
    }

    /// The rescale arm is exact: scaling both arguments by `2^k` changes no bit, through the window's
    /// edges and the cold arm on either side, while neither scaled argument is subnormal.
    #[test]
    fn scaling_both_arguments_by_a_power_of_two_changes_no_bit() -> Result<(), String> {
        let mut draw = Draw(2);
        for _ in 0..2_000 {
            let (y, x) = (draw.wide(-4, 4), draw.wide(-4, 4));
            let (want64, want32) = (at64(y, x)?, at32(y as f32, x as f32)?);
            for k in (-1015..=1015).step_by(7) {
                let s = libm::ldexp(1.0, k);
                assert!(
                    same64(at64(y * s, x * s)?, want64),
                    "({y:e}, {x:e}) at 2^{k}"
                );
            }
            for k in -118..=120 {
                let s = libm::ldexpf(1.0, k);
                let (y, x) = (y as f32, x as f32);
                assert!(
                    same32(at32(y * s, x * s)?, want32),
                    "({y:e}, {x:e}) at 2^{k}"
                );
            }
        }
        Ok(())
    }

    /// Monotone in `y` at fixed `x`, a thousand ulps either side of every place two formulas
    /// meet: binary64's midpoints between centres, `n/d = (2j + 1)/16`, and the swap, on both sides
    /// of the diagonal; binary32's fold at `½` and swap.
    #[test]
    fn monotone_where_the_formulas_meet() -> Result<(), String> {
        let mut at: Vec<f64> = (0..8).map(|j| f64::from(2 * j + 1) / 16.0).collect();
        at.push(1.0);
        at.extend((0..8).map(|j| 16.0 / f64::from(2 * j + 1)));
        let mut draw = Draw(9);
        for x in [1.0, draw.wide(-40, 40).abs(), draw.wide(-40, 40).abs()] {
            for &t in &at {
                let mut y = t * x;
                for _ in 0..1000 {
                    y = f64::from_bits(y.to_bits() - 1);
                }
                let mut last = at64(y, x)?;
                for _ in 0..2000 {
                    y = f64::from_bits(y.to_bits() + 1);
                    let next = at64(y, x)?;
                    assert!(
                        last <= next,
                        "f64 at ({y:e}, {x:e}): {last:e} then {next:e}"
                    );
                    last = next;
                }
            }
        }
        for at in [0.5f32, 1.0, 2.0] {
            let mut y = at;
            for _ in 0..1000 {
                y = f32::from_bits(y.to_bits() - 1);
            }
            let mut last = at32(y, 1.0)?;
            for _ in 0..2000 {
                y = f32::from_bits(y.to_bits() + 1);
                let next = at32(y, 1.0)?;
                assert!(last <= next, "f32 at y = {y:e}: {last:e} then {next:e}");
                last = next;
            }
        }
        Ok(())
    }

    /// The binary32 kernel's error in its own ulps, against `libm`'s binary64 `atan2` (whose error is
    /// `2^-29` of a binary32 ulp), over every octant, both folds and arguments `2^±60` apart.
    fn f32_error(
        draws: usize,
        seed: u64,
        subject: impl Fn(f32, f32) -> Result<f32, String>,
    ) -> Result<(f64, f32, f32), String> {
        let mut draw = Draw(seed);
        let mut worst = (0.0, 0.0f32, 0.0f32);
        for i in 0..draws {
            let (y, x) = if i % 2 == 0 {
                (draw.wide(-60, 60), draw.wide(-60, 60))
            } else {
                (draw.wide(-2, 2), draw.wide(-2, 2))
            };
            let (y, x) = (y as f32, x as f32);
            let want = libm::atan2(f64::from(y), f64::from(x));
            let got = f64::from(subject(y, x)?);
            let ulp = libm::ldexp(1.0, libm::ilogb(want).max(-126) - 23);
            let err = (got - want).abs() / ulp;
            if err > worst.0 {
                worst = (err, y, x);
            }
        }
        Ok(worst)
    }

    #[test]
    #[ignore = "a measurement, printed for 0105 (draft) to cite"]
    #[allow(clippy::print_stdout)]
    fn the_binary32_kernel_reads_its_error() -> Result<(), String> {
        let (err, y, x) = f32_error(1 << 20, 3, at32)?;
        println!("f32 max {err:.4} ulp at atan2({y:e}, {x:e})");
        let (err, y, x) = f32_error(1 << 20, 3, |y, x| Ok(libm::atan2f(y, x)))?;
        println!("libm::atan2f max {err:.4} ulp at atan2({y:e}, {x:e})");
        Ok(())
    }

    /// Binary64 against `libm`: a smoke bound, since `libm` is itself up to an ulp off; the measured
    /// readings are the corpus's (`helicoid:portable-atan2` over `real_atan2`).
    #[test]
    fn binary64_is_within_two_ulps_of_libm() -> Result<(), String> {
        let mut draw = Draw(4);
        let mut worst = 0u64;
        for _ in 0..1_000_000 {
            let (y, x) = (draw.wide(-30, 30), draw.wide(-30, 30));
            let (got, want) = (at64(y, x)?, libm::atan2(y, x));
            worst = worst.max(got.to_bits().abs_diff(want.to_bits()));
        }
        assert!(worst <= 2, "{worst} ulps");
        Ok(())
    }

    /// The guard every run keeps: within `0.51` ulp of the reference at binary64 and `0.5001` at
    /// binary32, near correct rounding, on a sample the ignored readouts take sixteen times over.
    #[test]
    fn the_kernels_stay_near_correct_rounding() -> Result<(), String> {
        let (err, y, x) = f64_error(1 << 16, 7, (-30, 30), at64)?;
        assert!(err < 0.51, "f64 {err} ulp at atan2({y:e}, {x:e})");
        let (err, y, x) = f32_error(1 << 16, 7, at32)?;
        assert!(err < 0.5001, "f32 {err} ulp at atan2({y:e}, {x:e})");
        Ok(())
    }
}
