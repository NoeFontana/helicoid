# The single task surface for helicoid. CI invokes these recipes 1:1.
# Recipes for the instrument (corpus, conformance, thresholds, oracles, envelope, determinism,
# bench) land with the `docs/PHASE1.md` step that owns each; none is stubbed to pass.

default:
    @just --list

# Build the whole workspace.
build:
    cargo build --workspace --all-targets

# Unit + integration tests (nextest) in the dev profile, then in release, where `debug_assert!` is
# compiled out and the out-of-domain behaviour the docs state is what is tested; then doctests.
test:
    cargo nextest run --workspace --no-tests=pass
    cargo nextest run --workspace --release --no-tests=pass
    cargo test --doc --workspace
    cargo nextest run -p helicoid-linalg --features mint
    cargo test --doc -p helicoid-linalg --features mint

# fmt, clippy with warnings denied, then the repository checks of `cargo xtask lint`. The
# workspace pass unifies `xtask`'s hidden `__sweep` into `helicoid`, so `helicoid` is linted again at
# its own default features, as a consumer builds it: an item only the sweep or the tests use is dead
# code there.
lint:
    cargo fmt --all -- --check
    cargo clippy --workspace --all-targets -- -D warnings
    cargo clippy -p helicoid -- -D warnings
    cargo xtask lint
    cargo clippy -p helicoid-linalg --all-targets --features mint -- -D warnings

# cargo-deny: advisories, licenses, sources for the whole workspace; the `0007` bans for the normal
# dependencies of everything but `xtask` (dev-dependencies such as `proptest` are out of scope).
audit:
    cargo deny check advisories licenses sources
    cargo deny --exclude-dev --exclude xtask check bans

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
    cargo "+$want" build -p helicoid-linalg --features mint --lib --locked

# No `std`, no `alloc`: the library crates on a bare-metal target. `--features __sweep` builds the
# coefficient kernel, which only it and the tests build until SO(3) uses it.
no-std:
    cargo build --target thumbv7em-none-eabihf -p helicoid-linalg -p helicoid --locked
    cargo build --target thumbv7em-none-eabihf -p helicoid --features __sweep --locked
    cargo build --target thumbv7em-none-eabihf -p helicoid-linalg --features mint --locked

# Library crates on wasm32, the coefficient kernel included as under `no-std`. Runs under wasmtime
# once the conformance subject exists.
wasm:
    cargo build --target wasm32-wasip1 -p helicoid-linalg -p helicoid --locked
    cargo build --target wasm32-wasip1 -p helicoid --features __sweep --locked
    cargo build --target wasm32-wasip1 -p helicoid-linalg --features mint --locked

# In-process subjects over the corpus: forward error per (fn, stratum) into conformance/results/,
# the table by max_u descending; fails on any non-finite output and when nothing was scored.
# Arguments: --subject NAME, --fn ID, --precision f64|f32 (f32 scores the `@f32` strata alone;
# default f64); or --self-test alone: the correct seeded kernel must fire no mechanism and every
# planted defect must fire its own (`docs/PHASE1.md` §10).
conformance *args:
    cargo xtask conformance {{args}}

# The threshold sweep at `f64` and `f32` (`docs/PHASE1.md` §6), of two targets: the seeded kernels
# (conformance/sweeps/thresholds-seeded.csv, xtask/src/seeded/generated.rs) and the arms `helicoid`
# ships, through its hidden `__sweep` feature (conformance/sweeps/thresholds.csv,
# crates/helicoid/src/coeffs/generated.rs); `cargo xtask thresholds [--check] [seeded|helicoid]`
# names one. Each generated file is compiled into xtask: a hand edit that no longer compiles stops
# this recipe too, and `git restore` of that file is the way back. When the corpus's series change,
# the first run writes the new series beside placeholder switches and fails: run it again.
thresholds:
    cargo xtask thresholds

# A fresh sweep must equal every committed file byte for byte; nothing is written.
thresholds-check:
    cargo xtask thresholds --check

# Regenerate the mpmath corpus into conformance/corpus (`docs/PHASE1.md` §4; needs `uv`).
corpus:
    cd conformance/generate && uv run --frozen python -m gen all --out ../corpus

# The generator's unit tests (stdlib unittest).
corpus-test:
    cd conformance/generate && uv run --frozen python -m unittest discover -s tests

# Regenerate into a temporary directory and compare with the committed corpus byte for byte:
# any differing, extra or missing file fails. Then the generator's own tests.
corpus-check: && corpus-test
    #!/usr/bin/env bash
    set -euo pipefail
    tmp=$(mktemp -d)
    trap 'rm -rf "$tmp"' EXIT
    (cd conformance/generate && uv run --frozen python -m gen all --out "$tmp")
    diff -rq "$tmp" conformance/corpus
