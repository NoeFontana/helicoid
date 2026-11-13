# helicoid-spline

Uniform cumulative B-splines of order `K = 2..6` on any `helicoid` Lie group: value, body velocity
and body acceleration by Sommer et al.'s recursions, with `O(K)` Jacobians on every control, either
side. `no_std`, no `alloc`, no `unsafe`.

The caller hands over a window of `K` controls, the segment parameter `u` and `1/h`; the crate owns
no storage and no time. Formulas: `docs/NUMERICS.md` §21. Decisions: `docs/decisions/0112`–`0116`.
