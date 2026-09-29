//! The shape of a generated switch (`docs/PHASE1.md` §6), defined next to `generated.rs`, which
//! `cargo xtask thresholds` writes with one per coefficient. Phase 3 defines its own beside
//! `helicoid::coeffs`; this one is the seeded kernels'.

/// `M` series terms below the branch variable `below`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Switch<S, const M: usize> {
    /// The switch in `z = θ²`: the series arm is taken where `z < below`, the exact arm elsewhere.
    pub(crate) below: S,
    /// The series in `z`, lowest power first.
    pub(crate) series: [S; M],
}

impl<S: Copy, const M: usize> Switch<S, M> {
    /// The switch and the series as a slice, whatever `M` is.
    pub(crate) fn parts(&'static self) -> (S, &'static [S]) {
        (self.below, &self.series)
    }
}
