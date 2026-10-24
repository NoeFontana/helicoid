# The single task surface for helicoid. CI invokes these recipes 1:1.
# Recipes for the instrument (corpus, conformance, thresholds, oracles, envelope, determinism,
# bench) land with the `docs/PHASE1.md` step that owns each; none is stubbed to pass.

default:
    @just --list

# Build the whole workspace.
build:
    cargo build --workspace --all-targets

# Unit + integration tests (nextest) in the dev profile, then in `ndebug` (`Cargo.toml`), where
# `debug_assert!` is compiled out and the out-of-domain behaviour the docs state is what is tested;
# then doctests.
test:
    cargo nextest run --workspace --no-tests=pass
    cargo nextest run --workspace --cargo-profile ndebug --no-tests=pass
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
    # The `ndebug` profile too. `debug_assertions` is off there, so every `cfg(not(debug_assertions))`
    # test body -- the ones `just test`'s `ndebug` run is the only thing that executes -- is cfg'd out
    # of the dev pass above and would otherwise never be checked against the denied `unwrap_used`,
    # `expect_used`, `panic`, `todo`, `unimplemented` and `dbg_macro`. The `mint` pass below stays
    # dev-only, as `just test`'s `mint` run is.
    cargo clippy --workspace --all-targets --profile ndebug -- -D warnings
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
# The sophus-rs oracle (`docs/PHASE1.md` §7): the runner's own checks, since it is
# workspace-excluded and `lint`, `test`, `doc` and `audit` never reach it, then the corpus through
# it, scored with the harness's exact metric. A runner that fails, or stops answering an id it
# owes, fails the recipe; a wrong or non-finite answer is recorded in its row, not a failure
# (oracles may be wrong: §7). `subject_version` is the pinned `sophus_lie`. Audited with the
# runner's own `deny.toml`, whose one advisory exception is the crate's, not the workspace's:
# cargo-deny finds it beside the manifest, so no `--config`, which cargo-deny 0.20 moved out of
# `check`. Rows in conformance/results/sophus_rs.csv.
oracle-sophus-rs:
    cargo fmt --manifest-path runners/sophus_rs/Cargo.toml -- --check
    cargo clippy --manifest-path runners/sophus_rs/Cargo.toml --locked --all-targets -- -D warnings
    RUSTDOCFLAGS='-D warnings' cargo doc --manifest-path runners/sophus_rs/Cargo.toml --locked --no-deps --document-private-items
    cargo deny --manifest-path runners/sophus_rs/Cargo.toml check advisories licenses sources
    cargo nextest run --manifest-path runners/sophus_rs/Cargo.toml --locked
    cargo xtask conformance --oracle sophus_rs

# The nalgebra oracle (`docs/decisions/0056` decision 4): as `oracle-sophus-rs`, for
# `runners/nalgebra` (`subject_version` is the pinned `nalgebra`), audited with the runner's own
# `deny.toml`, which cargo-deny finds beside the manifest: the root's bans `nalgebra`. It answers
# `eig3`, `chol_*`, `chol_solve_*`, `solve_cubic`, `quat_renormalize`, `mat2_inverse_adj`, `svd3`
# and `nearest_rotation` (`docs/decisions/0101`, `0102`), in binary64 and, on their `@f32`
# strata, in binary32 (`docs/decisions/0077`). Rows in conformance/results/nalgebra.csv and nalgebra--f32.csv.
oracle-nalgebra:
    cargo fmt --manifest-path runners/nalgebra/Cargo.toml -- --check
    cargo clippy --manifest-path runners/nalgebra/Cargo.toml --locked --all-targets -- -D warnings
    RUSTDOCFLAGS='-D warnings' cargo doc --manifest-path runners/nalgebra/Cargo.toml --locked --no-deps --document-private-items
    cargo deny --manifest-path runners/nalgebra/Cargo.toml check advisories licenses sources
    cargo nextest run --manifest-path runners/nalgebra/Cargo.toml --locked
    cargo xtask conformance --oracle nalgebra
    cargo xtask conformance --oracle nalgebra --precision f32

# The Sophus C++ oracle (`docs/decisions/0072`): the image of docker/oracles/sophus, pinned by base
# digest, Ubuntu snapshot and Sophus commit. Its build compiles the runner with warnings denied and
# runs the conversions' hand-computed tests; the image then runs offline over the corpus, as the
# owner of the answers' directory. Needs docker. `subject_version` names Sophus, Eigen, gcc and
# glibc; it answers in binary64, and `sim3_*`'s `@f32` strata with `Sophus::Sim3<float>`
# (`docs/decisions/0077`), `so3_from_matrix`'s with `Sophus::SO3<float>` (`0085`). `eig3` is the
# image's Eigen, `SelfAdjointEigenSolver<Matrix3<T>>`, at both precisions (`docs/decisions/0083`),
# `svd3` its `JacobiSVD<Matrix3<T>>` (`docs/decisions/0101`) and `nearest_rotation`
# `Sophus::SO3<T>::fitToSO3` (`docs/decisions/0102`). Rows in conformance/results/sophus.csv
# and sophus--f32.csv.
oracle-sophus:
    cargo xtask conformance --oracle sophus
    cargo xtask conformance --oracle sophus --precision f32

# The GTSAM C++ oracle (`docs/decisions/0074`): the image of docker/oracles/gtsam, pinned by base
# digest, Ubuntu snapshot and GTSAM commit, the library compiled in the image under the runner's
# flags (the first build takes minutes). Its build runs the conversions' hand-computed tests; the
# image then runs offline over the corpus, as the owner of the answers' directory. Needs docker.
# `subject_version` names GTSAM, Eigen, gcc and glibc; it answers in binary64. Rows in
# conformance/results/gtsam.csv.
oracle-gtsam:
    cargo xtask conformance --oracle gtsam

# The manif oracle (`docs/decisions/0075`): as `oracle-sophus`, for docker/oracles/manif, pinned by
# base digest, Ubuntu snapshot and manif commit; its build runs the conversions' hand-computed
# tests. `subject_version` names manif, Eigen, gcc and glibc; it answers in binary64, and in
# binary32 the one id with `@f32` strata it answers, `se3_inv_act` (`0103`). Rows in
# conformance/results/manif.csv and manif--f32.csv.
oracle-manif:
    cargo xtask conformance --oracle manif
    cargo xtask conformance --oracle manif --precision f32

# Every oracle runner (`docs/PHASE1.md` §3, §7): the two crates' checks and rows, then the three
# containers', each image built from its pins. Local only: it needs docker, and the first GTSAM
# build takes minutes.
oracles: oracle-sophus-rs oracle-nalgebra oracle-sophus oracle-gtsam oracle-manif

# The deep tier (`docs/decisions/0088`): every stratum whose sample count is below 64 and whose draw
# extends without moving its first samples, raised to 64, into the git-ignored
# conformance/corpus-deep/. Each committed stratum is checked to be its deep namesake's prefix; the
# manifest is committed as conformance/MANIFEST-deep.json. Local only, like the oracles and the
# envelope that read it.
corpus-deep:
    cd conformance/generate && HELICOID_CORPUS_DEEP=1 uv run --frozen python -m gen deep --committed ../corpus --out ../corpus-deep
    cp conformance/corpus-deep/MANIFEST.json conformance/MANIFEST-deep.json

# Draw the deep tier again into a temporary directory: its manifest must be the committed one, and
# a local conformance/corpus-deep/, where there is one, the same bytes.
corpus-deep-check:
    #!/usr/bin/env bash
    set -euo pipefail
    tmp=$(mktemp -d)
    trap 'rm -rf "$tmp"' EXIT
    (cd conformance/generate && HELICOID_CORPUS_DEEP=1 uv run --frozen python -m gen deep --committed ../corpus --out "$tmp")
    diff -q "$tmp/MANIFEST.json" conformance/MANIFEST-deep.json
    if [ -d conformance/corpus-deep ]; then diff -rq "$tmp" conformance/corpus-deep; fi

# The candidate and every oracle over the deep tier, into conformance/results-deep/ (`0088`), for
# `just envelope --deep`, which judges there the strata the committed tier cannot resolve. The
# runners as `just oracles` built them; needs docker.
conformance-deep:
    cargo xtask conformance --deep
    cargo xtask conformance --deep --precision f32
    cargo xtask conformance --oracle sophus_rs --deep
    cargo xtask conformance --oracle nalgebra --deep
    cargo xtask conformance --oracle nalgebra --deep --precision f32
    cargo xtask conformance --oracle sophus --deep
    cargo xtask conformance --oracle sophus --deep --precision f32
    cargo xtask conformance --oracle gtsam --deep
    cargo xtask conformance --oracle manif --deep
    cargo xtask conformance --oracle manif --deep --precision f32

# The one-variable twins (`0037`, draft): a candidate's program at a scalar whose transcendentals
# are the host's `std` instead of the `libm` crate, into conformance/results/*-host-std.csv. Both
# are planted subjects, so a plain `conformance` run skips them; `envelope` reads the rows of the
# candidate's *own* twin — `helicoid:host-std` for the library, the seeded stand-in's for the
# stand-in — to say whether a domination failure is D16's cost or is not explained by D16 either.
# `helicoid:host-std` is the shipped program at that scalar, so it answers every id the candidate
# answers; `seeded:host-std` answers `PHASE1.md` §10's kernels.
conformance-twin:
    cargo xtask conformance --subject seeded:host-std
    cargo xtask conformance --subject helicoid:host-std

# A paired comparison of two in-process subjects (`0038` decision 4, `0078`): per (fn, stratum),
# per id and over every id, the records each is better on, an exact sign test, the geometric mean
# and extremes of the per-record error ratio, and the ratio of the two maxima. Printed only: no
# result CSV, no baseline, no bar reads it. Arguments: --a NAME --b NAME [--fn ID]
# [--precision f64|f32], e.g. `just compare --a helicoid:homogeneous-r --b helicoid`.
compare *args:
    cargo xtask compare {{args}}

# The bars of `docs/PHASE1.md` §8 over conformance/results: domination over the best oracle and exact
# no-regress against conformance/baseline, per (fn, stratum, precision) on the max, then coverage of
# the corpus. Domination fails a maximum above the best oracle's only where the paired swap test on
# the records resolves it at p <= 0.001, and lists the others with their p; a stratum under 10
# records is unresolvable, and `--deep` judges it on `just conformance-deep`'s rows
# (`docs/decisions/0088`). It first reruns the in-process subjects at both precisions, so the candidate's rows are
# current (`docs/decisions/0077`: an oracle's `f32` rows pair with the candidate's `@f32` ones); the
# oracles' rows are what `oracle-sophus-rs`, `oracle-nalgebra`, `oracle-sophus`, `oracle-gtsam`, `oracle-manif` last wrote. `--check` also fails on
# a baseline or docs/evidence/ENVELOPE.md (its oracle columns aside) that `--bless` would write
# differently (an improved max is blessed in its own PR); `--bless` writes both, and nothing while
# domination, coverage or a dropped stratum fails (`--dry-run`: says what it would write). While the
# `helicoid` subject is not registered (Phase 3) it has no rows to judge: the run says so and checks
# coverage only; once it is, no rows fails. `--candidate seeded:correct` reads the seeded kernel as
# a stand-in, read-only in practice. It also reruns `conformance-twin`, so the attribution of a
# domination failure is never read from a stale twin.
envelope *args: conformance conformance-twin
    cargo xtask conformance --precision f32
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

# D16 (`docs/PHASE6.md` §3): the shipped subject's outputs over the whole corpus, one SHA-256 per
# id and precision, digested on this host and under wasmtime for `wasm32-wasip1`, then compared; any
# difference fails and names the id. The host digests twice, at `release` and at `ndebug`'s
# opt-level 0: `just test` runs at the second, so an optimizer that moved a bit would fail here.
# CI digests x86_64 (both profiles), aarch64 and wasm32 in separate jobs with the two recipes below
# and compares the four files with `cargo xtask determinism --compare`.
determinism: (determinism-digest "target/determinism/native.txt") (determinism-digest "target/determinism/ndebug.txt" "ndebug") (determinism-wasm "target/determinism/wasm32.txt")
    cargo xtask determinism --compare target/determinism/native.txt target/determinism/ndebug.txt target/determinism/wasm32.txt

# This host's digest table into `out`, built at `profile`.
determinism-digest out profile="release":
    mkdir -p "$(dirname '{{out}}')"
    cargo run --quiet --profile '{{profile}}' --locked -p xtask -- determinism --write '{{out}}'

# The `wasm32-wasip1` digest table into `out` (needs `wasmtime`). The guest reads the corpus through
# the checkout's absolute path, which is the one the harness was compiled with.
determinism-wasm out:
    mkdir -p "$(dirname '{{out}}')"
    cargo build --quiet --release --locked --target wasm32-wasip1 -p xtask
    wasmtime run --dir '{{justfile_directory()}}' target/wasm32-wasip1/release/xtask.wasm \
        determinism --write '{{justfile_directory()}}/{{out}}'

# The threshold sweep of the seeded kernels over the corpus (`docs/PHASE1.md` §6): writes
# conformance/sweeps/thresholds.csv and xtask/src/seeded/generated.rs. `coeffs/generated.rs`
# joins them with Phase 3. The second is compiled into xtask: a hand edit that no longer compiles
# stops this recipe too, and `git restore xtask/src/seeded/generated.rs` is the way back.
thresholds:
    cargo xtask thresholds

# A fresh sweep must equal every committed file byte for byte; nothing is written.
thresholds-check:
    cargo xtask thresholds --check

# The portable `atan2` candidate's tuning (`0105`, draft): binary64's centres and terms of P,
# binary32's terms of P. The kernels read the shape back from the fit at compile time.
atan_fit_tuning := "--f64 8:5 --f32 8"

# The portable `atan2` candidate's polynomials and tables (`0105`, draft) into
# conformance/fits/atan2.json, read by `helicoid:portable-atan2`.
atan-fit:
    cd conformance/generate && uv run --frozen python -m fit atan2 {{atan_fit_tuning}} --out ../fits/atan2.json

# Check conformance/fits/atan2.json is what `atan-fit` writes.
atan-fit-check:
    #!/usr/bin/env bash
    set -euo pipefail
    tmp=$(mktemp)
    trap 'rm -f "$tmp"' EXIT
    (cd conformance/generate && uv run --frozen python -m fit atan2 {{atan_fit_tuning}} --out "$tmp" > /dev/null)
    diff -u conformance/fits/atan2.json "$tmp"

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
