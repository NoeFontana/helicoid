# helicoid — Phase 7 Requirements: Continuous Time (`helicoid-spline`)

> **GATED by [`0011`](./decisions/0011-continuous-time-waits-for-a-consumer.md), not scheduled.**
> This is a requirements artifact, not an implementation authorization. §3 is questions, not
> answers; a question answered from this document instead of from a consumer's record is the
> failure 0011 prevents.

## 0.0 Gate

| Gate | Status |
|---|---|
| A consumer record exists naming the problem (one of the triggers below) | **Unmet** |
| That record states: group(s), knot spacing, order, derivative orders, precision, problem size, latency budget | **Unmet** |
| That record answers §3.2 (Kalibr-parity parameterization: needed or not) | **Unmet** |
| The core owes nothing new (§1 all delivered) | Met on Phase 3 completion |

Triggers, any one of: **locus-calib** camera–IMU calibration (continuous-time trajectory);
**locus_fusion** continuous-time or rolling-shutter estimation; **`tf_tree` Phase 6** (cumulative
B-spline interpolation, `tf_tree` 0009 §4).

## 1. What the core already owes, and delivers regardless

`Ad`, `ad`, `Jr`, `Jr⁻¹`, `Jl`, `Jl⁻¹` on every group including products (Phase 3); geodesics and
their Jacobians (Phase 4); `Dual` through all of it (Phase 2). The Sommer et al. (CVPR 2020)
recursions need exactly these. No spline-shaped hook is added to the core in advance.

## 2. Requirements sketch (for the record that opens the gate)

- Crate `helicoid-spline`, same budget as `helicoid` (0007).
- Uniform cumulative B-splines of order `K` (const generic), cumulative basis matrix $\tilde M$
  computed at compile time.
- Generic over the group, so SE(3) and split `Product<SO3, R3>` are two instantiations, not two
  code paths.
- Value, body velocity and body acceleration by the Sommer et al. recursions; control-point
  Jacobians in $O(K)$ per evaluation; time Jacobians.
- Verification: `Dual` through the spline; corpus ids from mpmath evaluating the definition
  ($\prod \mathrm{Exp}(\tilde B_j(u)\,d_j)$); oracle runners for basalt-headers and sophus-rs splines;
  benches beside Basalt.

## 3. Questions, not answers

1. **Non-uniform knots?** Kalibr and most VIO use uniform; rolling-shutter work sometimes does not.
2. **Kalibr parity.** Kalibr's rotation is a Euclidean spline on a rotation-vector
   parameterization, not a Lie-group spline. If locus-calib's baseline comparison needs parity,
   is that a `Chart`-parameterized Euclidean spline here, or does locus-calib own the confound in
   its evaluation protocol?
3. **`tf_tree`'s arena spline region** (`tf_tree` `PROJECT.md` §5.1): does a spline need stored
   control points in the arena, or does `tf_tree` evaluate from stored samples?
4. **`f32`?**
5. **Default for SE(3): cumulative on SE(3) or split?** The literature (Ovrén & Forssén; Sommer et
   al.) and Basalt's choice must be re-checked against the consumer's metric, not assumed.
6. **Gaussian-process priors (STEAM-style) instead of, or beside, splines?**
