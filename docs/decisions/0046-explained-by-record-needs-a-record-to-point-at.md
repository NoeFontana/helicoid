# 0046: "Explained by record" needs a record to point at

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** none yet; unblocks `PHASE3.md` §10

## Context

Two NORMATIVE passages give the domination bar an escape hatch, in nearly the same words:

- `PHASE3.md` §10: "**every paired stratum must be dominated**. A stratum where an oracle wins
  blocks the phase: fix the implementation **or write a record explaining why the oracle is
  measuring something else**. Bless the baseline; commit `docs/evidence/ENVELOPE.md`."
- `PHASE4.md` §5.2: "A stratum where `tf_tree_math` wins blocks the migration **until fixed or
  explained by record**."

**Neither hatch exists.** `cargo xtask envelope --bless` writes nothing while any domination
failure stands, and there is nowhere to record an explanation, so the two specs describe a state
the tool cannot represent. The consequence is not theoretical: `PHASE3.md` §0.0's envelope row
reads "Not started … fails on no baseline and 99 domination failures", `docs/evidence/` is empty,
and **no baseline has ever been blessed**, so the 1529 scored rows carry no no-regress bar at all.
The strictest bar in the project is currently buying nothing, because a bar that cannot be
satisfied is not enforced — it is skipped.

Some of those rows are, on the evidence, not defects. Measured over the 164 `f64` strata paired
with oracle #1, `helicoid` is worse on 13, and three of the thirteen agree with `tf_tree_math` to
**ten significant digits**:

| id | stratum | `helicoid` | `tf_tree_math` |
|---|---|---|---|
| `so3_exp` | `theta:1e-5` | 0.4952463054990454 | 0.4952463054734329 |
| `sen3_log_n1` | `theta:1e-4` | 0.3967683750061985 | 0.3967652075302560 |

Both programs are under half an ulp of correctly rounded, and the margin is $5\times10^{-11}$
relative. [`0038`](./0038-a-program-comparison-is-not-a-bar.md) considered this and **declined to
loosen the bar** — rightly: its decision 1 keeps `0006`'s maximum as the gate, because a bar that
passes a program because it is typically fine is the failure mode `0006` was written against. So
the resolution cannot be a tolerance on the bar. It has to be the hatch the specs already promise.

## Decision

1. **A committed exception table**, `conformance/baseline/exceptions.toml`, one row per excepted
   `(fn, stratum, precision)`, each carrying the record that explains it and a one-line reason:

   ```toml
   [[exception]]
   fn = "so3_exp"
   stratum = "theta:1e-5"
   precision = "f64"
   record = "0046"
   reason = "ten-digit tie: both programs under 1/2 u, margin 5e-11 relative"
   ```

2. **It excepts domination and nothing else.** No-regress, a non-finite candidate output, a
   candidate row with nothing scored, an oracle over another record count and coverage are all
   unaffected and remain unconditional failures. The exception does not change the candidate's
   scored `max_u`, which is still written to the baseline and still has to not regress — so an
   excepted stratum is *watched more closely* than a dominated one, not less.
3. **An exception is a citation, and the lint checks it.** `cargo xtask lint` fails when a row's
   `record` does not exist or is not `ready` — reusing the status reader
   `lint/drafts.rs` already has. A `draft` record cannot except a bar, by
   [`0040`](./0040-a-draft-is-not-a-parking-space.md) ("a draft authorises nothing").
4. **A stale exception fails the run.** If an excepted stratum is dominated after all, `--bless`
   and `--check` both fail and name the row to delete. An exception is a debt with a test attached;
   it cannot outlive the defect it describes, which is the failure mode of every waiver list ever
   written.
5. **`docs/evidence/ENVELOPE.md` prints the exceptions as a section of its own** — the stratum, the
   two maxima, the margin, and the record — so the page a reader trusts says what it is not
   claiming. `--check` compares it byte for byte as it already does.
6. **This record is itself the explanation for the ten-digit ties.** They are excepted on the
   stated ground: both programs' maxima are below $\tfrac12\,u$ and the margin is $5\times10^{-11}$
   relative, which is smaller than the difference between the two programs' `libm` builds would be
   on another host, so the ordering is not a property of either kernel. `0038` decision 2's paired
   instrument is the right way to compare the two programs on these strata, and it reports a
   direction, not a winner. **Only those strata**: the exception table is not where the other ten
   paired losses go, and this record excepts none of them.

## Rationale

The three candidates were a tolerance, a rewording, and a mechanism.

**A tolerance on the bar** — "pass when within $\epsilon$ of the oracle" — is what `0038` rejected
and it is rejectable on its own terms: the number would be typed, which `0004`'s spirit resists; it
would apply to every stratum, including ones where a 1 % margin is a real defect; and a program
that is 1 % worse everywhere would pass silently, which is exactly the regression the maximum bar
exists to catch.

**Rewording §5.2 and §10** to drop the hatch would make the bar honest but unsatisfiable: with one
unexplained stratum, no baseline can ever be blessed, so the no-regress bar — the *other* half of
D8, and the one that catches tomorrow's mistake rather than today's — stays switched off forever.
That is a worse trade than any single stratum's verdict.

**The mechanism** costs one committed file and two checks, keeps the max bar exactly as `0006` and
`0038` want it, and makes the escape hatch auditable: every exception has a name, a number, a
reason and a test that deletes it. It also turns the current 99-failure wall into a list of
decisions somebody has to write down, which is the state the specs assumed all along.

## Consequences

- `PHASE3.md` §10's "write a record explaining why" and `PHASE4.md` §5.2's "or explained by record"
  become operations, and both sections gain a sentence naming the file.
- The first blessed baseline becomes reachable, so D8's no-regress bar switches on for 1529 rows.
- Every exception is visible in three places: the table, the lint, and the evidence page. None of
  them can drift from the others without failing a run.
- An exception is cheap to write and impossible to forget, which is the opposite of the usual
  waiver file. The cost is that `--bless` now has a second input to review in a PR; the evidence
  page's new section is what makes that review one glance.

## Implementation plan

1. This record, `PHASE3.md` §10's and `PHASE4.md` §5.2's sentences — verified by `just lint`.
2. The table's reader and the two bars' changes in `xtask/src/envelope/` — verified by the
   existing synthetic-result tests extended with: an excepted failure that blesses, a stale
   exception that fails, an exception citing a missing record, an exception citing a `draft` record,
   and an exception on a no-regress failure, which must **not** pass.
3. The lint check for the cited record's existence and status — verified by `lint`'s own
   `every_check_runs`, which asserts an exact violation count and so needs its planted case.
4. The evidence page's section — verified by `--bless` twice byte for byte and `--check` on a hand
   edit, as the page's existing tests do.
5. The three ten-digit rows excepted, citing this record — verified by `just envelope` reaching
   coverage with those three no longer failing and the other ten still failing.

## Open questions

None.

## Further work

1. The other ten paired losses to oracle #1 — `sen3_log_n1` at $2.25\times$ and $2.09\times$,
   `sen3_exp_n1` near π at $1.39\times$ down to $1.01\times$, `so3_log` at $1.11\times$ and
   $1.10\times$ — are fixes or their own records, not entries here.
   [`0039`](./0039-the-sweeps-grid-stops-below-its-own-optimum.md) addresses the switch-bound ones.
2. Whether an exception should carry an expiry (a version, or a phase) so the list is pruned on a
   schedule rather than only when a defect is fixed.
3. `0038` decision 2's paired instrument would give each excepted stratum a direction and a sign
   test, which is strictly more informative than the two maxima this record's `reason` quotes.
