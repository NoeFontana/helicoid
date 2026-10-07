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
    # The release profile too. `debug_assertions` is off there, so every `cfg(not(debug_assertions))`
    # test body -- the ones `just test`'s release run is the only thing that executes -- is cfg'd out
    # of the dev pass above and would otherwise never be checked against the denied `unwrap_used`,
    # `expect_used`, `panic`, `todo`, `unimplemented` and `dbg_macro`. The `mint` pass below stays
    # dev-only, as `just test`'s `mint` run is.
    cargo clippy --workspace --all-targets --release -- -D warnings
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

# Library crates on wasm32, the coefficient kernel included as under `no-std`. Build only: the
# in-process conformance subject exists, its run under wasmtime is owed (PHASE1 §3, Phase 6).
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
# The tf_tree_math oracle (`docs/PHASE1.md` §7, `0010`): the runner's own checks, since it is
# workspace-excluded and `lint`, `test`, `doc` and `audit` never reach it, then the corpus through it,
# scored with the harness's exact metric into conformance/results/tf_tree_math.csv
# (`subject_version` is the pinned tf_tree commit). A runner that fails, or stops answering an id
# it owes, fails the recipe; a wrong or non-finite answer is recorded in its row, not a failure
# (oracles may be wrong: §7). The first build fetches the pinned commit from GitHub.
oracle-tf-tree-math:
    cargo fmt --manifest-path runners/tf_tree_math/Cargo.toml -- --check
    cargo clippy --manifest-path runners/tf_tree_math/Cargo.toml --locked --all-targets -- -D warnings
    RUSTDOCFLAGS='-D warnings' cargo doc --manifest-path runners/tf_tree_math/Cargo.toml --locked --no-deps --document-private-items
    cargo deny --manifest-path runners/tf_tree_math/Cargo.toml check advisories licenses sources
    cargo nextest run --manifest-path runners/tf_tree_math/Cargo.toml --locked
    cargo xtask conformance --oracle tf_tree_math

# The sophus-rs oracle (`docs/PHASE1.md` §7): as `oracle-tf-tree-math`, for `runners/sophus_rs`
# (`subject_version` is the pinned `sophus_lie`), audited with the runner's own `deny.toml`, whose
# one advisory exception is the crate's, not the workspace's: cargo-deny finds it beside the
# manifest, so no `--config`, which cargo-deny 0.20 moved out of `check`. Rows in
# conformance/results/sophus_rs.csv.
oracle-sophus-rs:
    cargo fmt --manifest-path runners/sophus_rs/Cargo.toml -- --check
    cargo clippy --manifest-path runners/sophus_rs/Cargo.toml --locked --all-targets -- -D warnings
    RUSTDOCFLAGS='-D warnings' cargo doc --manifest-path runners/sophus_rs/Cargo.toml --locked --no-deps --document-private-items
    cargo deny --manifest-path runners/sophus_rs/Cargo.toml check advisories licenses sources
    cargo nextest run --manifest-path runners/sophus_rs/Cargo.toml --locked
    cargo xtask conformance --oracle sophus_rs

# The one-variable twins (`0037`, draft): a candidate's program at a scalar whose transcendentals
# are the host's `std` instead of the `libm` crate, into conformance/results/*-host-std.csv. Both
# are planted subjects, so a plain `conformance` run skips them; `envelope` reads the rows of the
# candidate's *own* twin — `helicoid:host-std` for the library, the seeded stand-in's for the
# stand-in — to say whether a domination failure no `libm`-crate oracle wins is D16's cost or is
# not explained by D16 either. `helicoid:host-std` is the shipped program at that scalar, so it
# answers every id the candidate answers; `seeded:host-std` answers `PHASE1.md` §10's kernels.
conformance-twin:
    cargo xtask conformance --subject seeded:host-std
    cargo xtask conformance --subject helicoid:host-std

# The bars of `docs/PHASE1.md` §8 over conformance/results: domination over the best oracle and exact
# no-regress against conformance/baseline, per (fn, stratum, precision) on the max, then coverage of
# the corpus. It first reruns the in-process subjects, so the candidate's rows are current; the
# oracles' rows are what `oracle-tf-tree-math`, `oracle-sophus-rs` last wrote. `--check` also fails on
# a baseline or docs/evidence/ENVELOPE.md (its oracle columns aside) that `--bless` would write
# differently (an improved max is blessed in its own PR); `--bless` writes both, and nothing while
# domination, coverage or a dropped stratum fails (`--dry-run`: says what it would write). While the
# `helicoid` subject is not registered (Phase 3) it has no rows to judge: the run says so and checks
# coverage only; once it is, no rows fails. `--candidate seeded:correct` reads the seeded kernel as
# a stand-in, read-only in practice. It also reruns `conformance-twin`, so the attribution of a
# domination failure is never read from a stale twin.
envelope *args: conformance conformance-twin
    cargo xtask envelope {{args}}

# The bench gate (`docs/PHASE1.md` §9): criterion over the shipped coefficient kernel, then the
# paired bootstrap 95% CI of the ratio per benchmark, failed when the whole CI lies above `1 + δ`.
# `just bench --against <bench-binary>` is the gate: it runs each benchmark as baseline, candidate,
# baseline, so `δ` is an A/A control measured beside the comparison and not read from a file.
# `just bench --aa --bless` logs what this host's A/A noise has shown into baseline/HOST.md, which
# says whether the machine is quiet enough to bother; it is a log and not the allowance, because the
# same protocol measured 0.0161 one session and 0.6864 the next (`docs/decisions/0033`, draft).
# Pin a core (`taskset -c N just bench ...`) and run it on a quiet machine.
bench *args:
    cargo xtask bench-gate {{args}}

# Re-run the gate's decision rule over a recording, measuring nothing (`docs/decisions/0035`, draft).
# `just bench --against <binary> --record <dir>` writes one; this replays it in milliseconds, so the
# statistic can be changed and re-judged without paying for a quiet machine twice.
bench-replay dir:
    cargo xtask bench-gate --replay {{dir}}

# The threshold sweep of the seeded kernels over the corpus (`docs/PHASE1.md` §6): writes
# conformance/sweeps/thresholds.csv and xtask/src/seeded/generated.rs. `coeffs/generated.rs`
# joins them with Phase 3. The second is compiled into xtask: a hand edit that no longer compiles
# stops this recipe too, and `git restore xtask/src/seeded/generated.rs` is the way back.
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

# Everything `cargo publish` checks, without publishing (`0049` step 2). The tree must be clean:
# no `--allow-dirty`, because a release is cut from committed state.
#
# The order is not a preference. `cargo package` resolves a `path` dependency that carries a
# `version` from the index, so `helicoid` cannot be packaged *at all* -- not even with
# `--no-verify` -- until `helicoid-linalg` is on crates.io. Until then its only available check is
# the file list, which is what this prints. Publish `helicoid-linalg` first, wait for the index,
# then rerun this and publish `helicoid`.
publish-check:
    #!/usr/bin/env bash
    set -euo pipefail
    cargo publish -p helicoid-linalg --dry-run
    if cargo package -p helicoid --no-verify >/dev/null 2>&1; then
        cargo publish -p helicoid --dry-run
    else
        echo
        echo "helicoid: dry-run deferred, helicoid-linalg is not on the index yet."
        echo "What would be packaged:"
        cargo package -p helicoid --list
    fi

# The exhaustive guards: every binary32 bit pattern, for the claims that license re-spelling a
# `sin_cos` as a `sin` or a `cos` (`0052`, `0053`). Minutes, so **not** in `just test` — two shard
# sets of eight, run concurrently, about four minutes of CPU each.
#
# Run it after `cargo update` touches `libm`: it is a caret dependency, the two functions of each
# pair are *not* the same expression, and `solve_cubic`'s committed golden bits rest on the identity.
# `just test` runs the dense samples, which would catch a structural change but not a rare one.
exhaustive:
    cargo nextest run --release -p helicoid-linalg \
        -E 'test(binary32_agrees_exhaustively) + test(cos_agrees_on_all_of_binary32)' \
        --run-ignored all
