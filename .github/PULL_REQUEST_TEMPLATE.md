<!-- Work scoped by a spec cites its section. Anything the specs do not cover starts as a
     `draft` record in docs/decisions/, not as a PR (CONTRIBUTING.md). -->

**Spec section / decision record:** <!-- e.g. docs/PHASE1.md §3, or 00NN -->

## What and why

## Checklist

- [ ] `just lint`, `just test`, `just audit` pass; `just no-std wasm msrv doc` if a library crate changed.
- [ ] No new dependency for a library crate beyond `0007`'s set.
- [ ] A formula, coefficient or switch point change is a `NUMERICS.md` edit **and** a record; generated files were regenerated, not edited.
- [ ] `CHANGELOG.md` `[Unreleased]` says what changed and whether it breaks anything.
- [ ] Each spec's §0.0 status table is updated if this changes what is implemented.
