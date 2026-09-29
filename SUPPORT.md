# Support policy

A **single-maintainer, pre-1.0 project**; the expectations below are the ceiling.

## Response expectations

| Kind of report | Expectation |
|---|---|
| Security vulnerability ([`SECURITY.md`](./SECURITY.md)) | Acknowledged within 7 days; fix or advisory within 90 |
| Wrong result above the oracle envelope, or a bit-identity break across targets | Triaged within 14 days |
| Panic or non-finite output on valid input | Triaged within 14 days |
| Any other bug | Best effort. No timeline promised |
| Feature request | Best effort, and likely declined; see *What is not supported* |
| Question | Best effort. `docs/` answers most |

"Triaged" means read and labelled, not fixed.

## What is supported

- **x86_64 and aarch64 Linux, and wasm32**, for the library crates (bit-identical outputs, D16).
- **The current release**, and only it: no backports.
- Nothing is published yet; the workspace is `0.0.x`, where every release may break every other
  ([`CHANGELOG.md`](./CHANGELOG.md)). Each spec's §0.0 status table says what is implemented.

## What is not supported

- **Anything a `docs/` status table marks not implemented.**
- **Feature requests that widen the scope** ([`0009`](./docs/decisions/0009-what-helicoid-does-not-own.md)):
  sensor models, solvers, preintegration, splines, SIMD lanes. Refute the recorded argument
  first, through [`docs/decisions/`](./docs/decisions/).

## MSRV policy

The minimum supported Rust version is **1.87**, declared in `[workspace.package] rust-version`
and never above the lowest MSRV among consumers (D17). `just msrv` builds at that version.
Suspended for the whole `0.0.x` line; the rule (a bump is a minor-version bump) returns at `0.1.0`.
**MSRV is raised only for a reason written in the raising commit.**

## Contributing

See [`CONTRIBUTING.md`](./CONTRIBUTING.md).
