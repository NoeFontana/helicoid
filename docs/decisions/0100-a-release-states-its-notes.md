# 0100: A release states its notes

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** this record's PR

## Context

`release.yml` publishes to crates.io on a `v*` tag and does nothing else.
- **No release page:** no GitHub release exists for `v0.0.1` to `v0.0.4`.
- **No check before publishing:** nothing tied the tag to the workspace's version, or the
  changelog to the tag.
- **What a reader cannot find:** what changed in a release, and which decision records it carries,
  lived only in `CHANGELOG.md` and the commit log.
- **The procedure was precedent:** it existed as the `release: 0.0.4` commit, not as text.

The owner asked for patch releases to carry automated GitHub release notes (2026-10-16).

## Decision

1. **`cargo xtask release-notes VERSION [--out FILE]`** writes the notes from what the repository
   already keeps:
   - `CHANGELOG.md`'s `## [VERSION]` section, verbatim, its relative links made absolute at the tag;
   - the decision records added or amended since the previous `v*` tag, with their titles and
     their `**Status:**` at the tag;
   - the crates, the MSRV and the exact pin.

   It fails where the tag is not `[workspace.package] version`, or the changelog has no non-empty
   section for it.
2. **`release.yml` runs three jobs in order.**
   - **`notes`** first, so a mis-tagged or unstamped release publishes nothing.
   - **`publish`**, as before.
   - **`github-release`**, which runs `gh release create --verify-tag --notes-file notes.md
     --generate-notes`. GitHub appends the merged pull requests, grouped by
     `.github/release.yml` (dependabot's under *Dependencies*).
3. **`CONTRIBUTING.md` *Releasing* states the procedure.**

## Rationale

- **One curated source.** The changelog is the half a human writes, and the record table and the
  pull-request list are derived. Nothing is written twice, so nothing drifts.
- **Fail before the irreversible step.** A crates.io version cannot be taken back. The notes job
  checks what a human gets wrong at release time, the tag and the stamp, before `publish` runs.
- **The record table carries the project's contract.** A `draft` record's change ships, but the
  record is not yet the contract, and the release page should say which.
- **Rejected: a hand-written GitHub release.** It is a second copy of the changelog, and the step
  most likely to be skipped.
- **Rejected: GitHub's generated notes alone.** They list pull-request titles; the changelog says
  what changed for a user and whether it breaks anything.

## Consequences

- **Every `v*` tag from `v0.0.5` gets a release page.** Earlier tags keep none.
- **A tag whose changelog section is missing fails the release,** and no crate is published.

## Implementation plan

1. `xtask/src/release_notes.rs` and its tests, including
   `the_committed_changelog_has_a_section_for_the_workspace_version`; `release.yml`;
   `.github/release.yml`; `CONTRIBUTING.md`. Verified by `cargo nextest run -p xtask`, `just lint`
   and the `v0.0.5` release itself.

## Open questions

None.
