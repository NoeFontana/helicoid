# 0040: A draft is not a parking space

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** the lifecycle in `decisions/README.md`, and the triage of four records below.

## Context

**Seventeen of thirty-eight records are `draft`, and they are the last twelve consecutively.**
Everything through [`0022`](./0022-real-owes-acos-and-cos.md) converged to `ready` or
`implemented`; of `0026` onward, only `0027` did. A draft authorises nothing
(`decisions/README.md`), so the project has spent four months unable to build on its own
measurements: [`0032`](./0032-domination-charges-helicoid-for-d16.md)'s 0.74 u, `0033`'s estimator
choices, `0037`'s resolution ceiling and `0039`'s 90× are all measured, reproducible, and
uncitable as settled.

The cause is not the evidence. It is the criterion, read more strictly than it is written.

- **`decisions/README.md` says `ready` means "every open question the *Decision* depends on is
  resolved".** In practice every open question has been treated as a blocker.
- **A record that records a measurement never runs out of open questions**, because a measurement's
  natural product is the next measurement. `0030` has nine, `0039` seven. Held to "no open
  questions", such a record can never leave `draft`, however settled its Decision.
- **[`0026`](./0026-write-dense-pays-for-zeros-not-for-checks.md) has zero open questions and a
  fully concrete Decision**, and has been `draft` regardless. Nothing was blocking it. That is the
  clearest evidence that the friction is process and not substance.
- **Two records are `draft` for a different reason and do not say so.**
  [`0028`](./0028-two-ready-specs-disagree-on-the-blocks-of-sen3jac.md) and
  [`0029`](./0029-a-product-needs-a-way-to-be-built-from-its-factors.md) both open with "Nothing is
  decided. A draft authorises nothing", state three options, and recommend one. They are not
  waiting for evidence; they are waiting for a choice. A reader cannot tell those two states apart
  from the status line, so neither gets made.

## Decision

1. **The criterion is applied as written.** A question the *Decision* does not depend on does not
   block `ready`. A record whose Decision is concrete and whose implementation plan is concrete is
   `ready` even if it ends with a list of things worth measuring next.
2. **Records separate the two kinds of question.** *Open questions* are the ones the Decision
   depends on, and they block; a new **Further work** section holds the rest, and it does not. A
   record with an empty *Open questions* section and a concrete Decision is `ready` by inspection,
   which is a property a reader can check in one line. `template.md` gains the section.
3. **A record that is waiting for a choice says so**, in its *Status* line: `draft (awaiting a
   decision)` rather than bare `draft`. "Awaiting evidence" and "awaiting a decision" are different
   states with different remedies, and conflating them is why `0028` and `0029` sat for months with
   their recommendation already written.
4. **A measurement is settled by being reproducible, not by being final.** Where a record's
   *Measurement* section states a protocol another session can rerun, that measurement is citable as
   measured once the record is `ready`, whatever the record proposes. This is what lets `0039` cite
   `0032`'s `atan2` column and `0038` cite `0039`'s arm timings without either being
   `implemented`.
5. **Triage applied now**, under the owner's instruction to decide, with each record's own
   recommendation accepted where it had one:

   | record | was | becomes | why |
   |---|---|---|---|
   | `0026` | draft | **ready** | zero open questions, concrete Decision; nothing was blocking it |
   | `0028` | draft | **ready**, option **A** | narrow `SEn3Jac` to `pub(crate)`, correct `PHASE3.md` §5 — the only reversible option, and no external consumer exists |
   | `0029` | draft | **ready**, option **A** | `from_parts`/`parts`, landing with the SO(3) PR — additive, and it agrees with `SE3::from_rt` one level down |
   | `0039` | draft | **ready** | its three limits are measured for accuracy *and* cost; its remaining questions are further work |

   The other thirteen drafts are not triaged here. Triaging a record one has not read is the failure
   mode this record exists to prevent.

## Rationale

The lifecycle was written for *decision* records, where "no open questions" is reachable, and the
instrument phase produced mostly *measurement* records, where it is not. `0006` is the reason —
"the instrument comes first" — so the drift was a consequence of following the roadmap, not of
ignoring it. The fix is to name the two kinds of question rather than to relax the bar: a Decision
that depends on an unresolved question must still stay `draft`, and item 1 does not change that.

Item 3 costs one parenthesis and buys the thing that was actually missing. A reader scanning
`grep -m1 -H '^\*\*Status:' docs/decisions/0*.md` — which `CLAUDE.md` tells an agent to run — sees
seventeen identical words today and cannot tell which of them is a question for them.

Item 4 is the one that compounds. A measurement that cannot be cited has to be redone, and this
session redid `0032`'s `atan2` swap as a subject precisely because the original script was not kept
— the right outcome for a tool, and pure waste if it had only needed a citation.

## Consequences

- **Four records leave `draft` and one spec row and one visibility change become implementable.**
  `0028` option A is a one-line narrowing plus a `PHASE3.md` §5 correction; `0029` option A waits
  for the SO(3) PR by its own decision; `0039`'s plan needs the sweep rerun, which is one command.
- `0013`'s "a **draft** and authorises nothing" sentence in `CLAUDE.md` stays true and stays
  necessary: item 1 widens what may be `ready`, not what a draft authorises.
- The `[drafts]` lint check keeps working unchanged — it reads the status line and requires the word
  `draft` beside a citation of one. Promoting a record makes its citations' "(draft)" stale rather
  than wrong, and the four promoted here are corrected in this change.
- **The remaining thirteen drafts are now a visible queue rather than an invisible one.** That is
  the point: the number was 17 and nobody had counted it.

## Implementation plan

1. `decisions/README.md`'s lifecycle gains items 1–4, `template.md` gains *Further work*, the four
   records' statuses change and their stale "(draft)" citations are corrected — verified by
   `cargo xtask lint` and by the status survey above. **This change.**
2. `0028` option A: `SEn3Jac`'s fields to `pub(crate)`, `PHASE3.md` §5's declaration corrected.
   **Owed, next.**
3. The other thirteen drafts triaged by their owner, each to `ready`, `draft (awaiting a
   decision)`, or superseded. **Owed.**

## Open questions

None. The Decision is about this repository's own process and depends on no measurement.

## Further work

- Should a *measurement* record have a kind of its own, rather than relying on item 4? Four of the
  thirteen remaining drafts record a measurement and propose nothing a reader would implement.
- The same finding has now arrived three times by three routes — `0033` for latency, `0038` for a
  program comparison, `0039` for the sweep's objective: **a scalar maximum discards the distribution
  it is a maximum of.** Three records state it locally and none states it as a principle. A
  successor should.
