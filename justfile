# The single task surface for helicoid. CI invokes these recipes 1:1.
# Recipes for the instrument (corpus, conformance, thresholds, oracles, envelope, determinism,
# bench) land with the `docs/PHASE1.md` step that owns each; none is stubbed to pass.

default:
    @just --list

# Build the whole workspace.
build:
    cargo build --workspace --all-targets

# Unit + integration tests (nextest), then doctests.
test:
    cargo nextest run --workspace --no-tests=pass
    cargo test --doc --workspace

# fmt, clippy with warnings denied, then the repository checks of `cargo xtask lint`.
lint:
    cargo fmt --all -- --check
    cargo clippy --workspace --all-targets -- -D warnings
    cargo xtask lint

# cargo-deny: advisories, licenses, sources. Library-crate bans (`0007`) are owed to `cargo xtask lint`.
audit:
    cargo deny check

# rustdoc with warnings denied.
doc:
    RUSTDOCFLAGS='-D warnings --cfg docsrs' cargo doc --no-deps --all-features --workspace

# The design book (`cargo install mdbook mdbook-katex`; NUMERICS.md is LaTeX-heavy).
book:
    mdbook build docs

# The declared floor, `--locked`. No fallback to stable: a missing toolchain stops the recipe.
msrv:
    #!/usr/bin/env bash
    set -euo pipefail
    want=$(grep -m1 '^rust-version' Cargo.toml | cut -d'"' -f2)
    test -n "$want" || { echo "no rust-version in Cargo.toml"; exit 1; }
    rustup toolchain list | grep -q "^$want" \
        || { echo "the floor is $want; install it: rustup toolchain install $want --profile minimal"; exit 1; }
    cargo "+$want" build --workspace --lib --bins --locked

# No `std`, no `alloc`: the library crates on a bare-metal target.
no-std:
    cargo build --target thumbv7em-none-eabihf -p helicoid-linalg -p helicoid --locked

# Library crates on wasm32. Runs under wasmtime once the conformance subject exists.
wasm:
    cargo build --target wasm32-wasip1 -p helicoid-linalg -p helicoid --locked
