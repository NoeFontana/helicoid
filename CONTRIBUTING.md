# Contributing

Agents also read [`CLAUDE.md`](./CLAUDE.md).

[`docs/PROJECT.md`](./docs/PROJECT.md) is the contract; read it, then
[`docs/NUMERICS.md`](./docs/NUMERICS.md) if you touch a formula, then the spec of the phase you
touch, **before** changing code. Each spec's §0.0 status table outranks every other statement.
When the documents do not answer a question, ask; do not invent a formula, coefficient, switch
point, tangent ordering or sign.

## Prerequisites

`rustup` (the toolchain and its targets come from `rust-toolchain.toml`; also `rustup toolchain
install 1.87 --profile minimal` for `just msrv`), [`just`](https://github.com/casey/just),
`cargo-nextest`, `cargo-deny`, and `mdbook` + `mdbook-katex` for `just book`. The reference
generator (`just corpus`, `just corpus-check`) needs [`uv`](https://docs.astral.sh/uv/).

## Workflow for significant changes

A change the specs do not cover (a public item, crate, dependency, convention, stratum family or
chart) starts as a **decision record**:

1. Copy [`docs/decisions/template.md`](./docs/decisions/template.md) to
   `docs/decisions/NNNN-kebab-case-title.md` (next number); status `draft`.
2. Open a PR with just the record; once its open questions are resolved, flip it to `ready`.
3. Implement under PRs that link the number, one per *Implementation plan* step.
4. When all merge, flip it to `implemented` and list the PR numbers.

A formula change is a `NUMERICS.md` edit **and** a record. Bug fixes, behaviour-preserving
refactors and dependency bumps need no record.

## Pull-request checklist

- [ ] `just lint test audit` pass; plus `just no-std wasm msrv doc` when a library crate changed.
- [ ] Library dependencies stay within [`0007`](./docs/decisions/0007-the-budget-a-foundation-can-afford.md).
- [ ] Generated files are regenerated, never hand-edited.
- [ ] Comments cite a symbol, never a line number, and state the decision and stop.
- [ ] `CHANGELOG.md` `[Unreleased]` says what changed and whether it breaks anything.
- [ ] The relevant §0.0 status table is updated.

Tests run under nextest: `cargo nextest run -p helicoid -- <name>`.

## License

Contributions are dual-licensed MIT OR Apache-2.0, as stated in [`NOTICE`](./NOTICE).
