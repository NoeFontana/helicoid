# 0045: Two Phase 4 checks cannot be taken as written, and the geodesic strata

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** plan steps 1 (this record and the `PHASE4.md` edits), 2 (`gen/geodesic.py`, `gen/check_geodesic.py`, the `geo:*` family, the registry entries and the generator's tests) and 3 (`just corpus`, the two committed files, `PHASE1.md` §4.3 and §4.4, the `coverage::OWED`-to-`required` move) have landed. Step 4, §3's invariance tests, waits for `LieGroup::geodesic` to exist

## Context

`PHASE4.md` §3 and §4 are both marked **NORMATIVE**, and neither can be executed as written. Both
findings come from `docs/maths/geodesics.md`, whose `Checked:` blocks give the protocol; the
`docs/maths/index.md` open-item list already carries both.

### §4's reference cannot be `mp.logm` near π

§4 defines `so3_geodesic` and `se3_geodesic` with "reference `mp.expm`/`mp.logm`" and a
`geo:near-pi` stratum at relative rotation $\pi - 10^{-k}$. Measured (SE(3), $X_0 = I$,
$t = \tfrac12$, 6 poses per $\theta$, 40 digits, **mpmath 1.3.0 and 1.4.1 identical**):
`mp.expm(t*mp.logm(Δ))` agrees with the geometric form to $9.2\times10^{-41}$ for
$\theta \le 3.0$, then returns a **complex, non-principal** logarithm — error $1.0$ to $1.9$ —
from $\theta = 3.03$ for 4 of 6 poses and for **all 6** at $\theta \ge 3.05$. So the stratum §4
requires is exactly the stratum its stated reference gets wrong, and by $O(1)$, not by ulps.

This is not new to the geodesic: the same defect is recorded against `so3_log`, `se2_log`,
`s2_local` and `sim3_log`, whose onsets fall with the extra parameter. It is a property of
`mp.logm` on a rotation, not of the formula.

### §3's right-invariance failure is vacuous under the law `Product` actually has

§3 requires, NORMATIVE: "**`Product<SO3, R3>` fails right-invariance, positively:** a fixed-seed
test asserts `max_err > 1e-6`. **Do not "fix" it.**" The reasoning is sound — the lerp is a
world-frame straight line while right-invariance forces the screw coupling
$\rho = J_l^{-1}(\varphi)\mathbf t$, and fixing it gives `ScLerp`, a different function. But
`Product` composes **componentwise**, $(R,t)(R',t') = (RR', t+t')$, and under its own law
GE.5(a) proves $\gamma$ is left- **and** right-invariant: measured at $7.3\times10^{-111}$
(rotation) and $1.7\times10^{-111}$ (translation) over 18 random pairs at 110 digits. A test
asserting `max_err > 1e-6` against `Product`'s own `Mul` cannot pass — it would be asserting a
falsehood, and the first green run would delete it as broken.

The failure is real only when $(R, t)$ is read as an SE(3) pose, $a\cdot H = (R_aR_H,\ R_at_H + t_a)$.
GE.5(b) gives it in closed form: the translation gap is $R_0M_s t_H$ with
$M_s = (1-s)I + sE - E_s$, of norm $\lvert\mu_s(\theta)\rvert\,\lVert t_{H\perp}\rVert$ where
$t_{H\perp}$ is the part of $t_H$ orthogonal to the axis and
$\lvert\mu_s\rvert = \tfrac12 s(1-s)\theta^2 + O(\theta^4)$. It vanishes **iff**
$s \in \{0,1\}$, $\theta = 0$, or $t_H$ is parallel to the axis — three ways a careless fixture
makes the NORMATIVE test vacuous a second time.

### A new stratum family is a record

`PHASE1.md` §4.4 and `gen/strata.py`'s own comment make a new stratum family a decision record
before it is code. §4 fixes the family's content; this record is where it is adopted.

## Decision

1. **`PHASE4.md` §4's reference becomes the geometric form, and `mp.logm` is named as excluded.**
   The reference is the geometric quaternion $\mathrm{Log}$ — `atan2` of the axis norm against $w$,
   as GE.2's own `Checked:` block computes it — with a $\mathsf V$/$\Gamma_1$ solve for the
   translation block, then `mp.expm` of $t\,d$. `mp.expm` stays; only `mp.logm` goes. The spec says
   why in one clause, so nobody reinstates it.
2. **The reference is cross-checked against `mp.logm` where `mp.logm` is sound** —
   $\theta \le 3.0$, where the two agree to $9.2\times10^{-41}$ — as the generation-time check
   (`gen/check_geodesic.py`), and the near-π cells are checked against the group identities instead
   (the endpoints $\gamma(X_0,X_1,0) = X_0$ and $\gamma(X_0,X_1,1) = X_1$, the symmetry
   $\gamma(X_0,X_1,t) = \gamma(X_1,X_0,1-t)$ of GE.2(c), and constant body velocity). A reference
   nobody can cross-check is a reference nobody should commit.
3. **The `geo:*` family is adopted as §4 states it** — `geo:consecutive`
   ($\lVert d\rVert \in [10^{-9}, 10^{-3}]$ with $\lVert\mathbf t_0\rVert$ up to $10^4$, the
   kilohertz edge), `geo:generic`, `geo:near-pi` ($\pi - 10^{-k}$), at
   $t \in \{0, 10^{-9}, 0.25, 0.5, 1-10^{-9}, 1\}$ plus uniform samples — with two additions the
   maths forces:
   - **`geo:near-pi` stops strictly below $\pi$ and records its margin.** GE.13(d): at
     $\theta = \pi$ there are two preimages, $\mathrm{Exp}\,d' = \mathrm{Exp}\,d$ to
     $3.9\times10^{-111}$, and the two geodesics through them differ at $t = \tfrac12$ by $0.83$ to
     $1.66$ — and $\gamma$ at $\theta = \pi \mp 10^{-9}$ differs by $1.659$, $0.829$, $1.354$ at
     $s = \tfrac14, \tfrac12, \tfrac34$. So two programs that disagree on the quaternion sign by
     one ulp of $\theta$ disagree on the answer by $O(1)$. The stratum is a measurement of
     conditioning, not of agreement, and the spec says which.
   - **$t = 0$ and $t = 1$ are kept and are exact rows.** GE.7(b) and GE.13(f): the exact arm at
     $s = 0$ is the identity with **zero** deviation, and at $s = 1$ within $2.0\,u$ per quaternion
     component and $4.9\,u$ in the translation relative to $\max(1,\lVert\mathbf t\rVert)$. These
     are the cheapest rows in the corpus and they catch an endpoint error no interior $t$ sees.
4. **`PHASE4.md` §3 says which law $a\cdot H$ uses.** The right-invariance row reads: under
   `Product`'s own componentwise law the geodesic **is** right-invariant (GE.5(a)) and the test is
   that it is; the positive failure is asserted against the **SE(3)** reading of $(R,\mathbf t)$,
   GE.5(b). The fixture is then constrained rather than fixed-seed-and-hope: $s = \tfrac12$,
   $t_H$ with a component orthogonal to the axis, and $\theta^2\lVert t_{H\perp}\rVert \gtrsim
   10^{-5}$ to clear `max_err > 1e-6` — the bound read off $\tfrac12 s(1-s)\theta^2\lVert
   t_{H\perp}\rVert$, and the test's comment carries that inequality so a later fixture edit cannot
   silently make it vacuous.
5. **The invariance bounds are measured per $\lVert\mathbf t_G\rVert$.** `geodesic(G·a, G·b, s)`
   forms $a^{-1}G^{-1}Gb$, which equals $a^{-1}b$ only to rounding, and §4's own strata put
   $\lVert\mathbf t_0\rVert$ at $10^4$. The left-invariance bound is recorded at
   $\lVert\mathbf t_G\rVert \in \{0, 1, 10^4\}$ and set from the measurement, per `0006`, not
   guessed from the $\lVert\mathbf t_G\rVert = 1$ case.

## Rationale

For §4 the alternative is to generate `geo:near-pi` from `mp.logm` anyway and treat the result as
the reference. It would produce a committed corpus whose near-π cells are wrong by $O(1)$, and the
first thing that happened would be a hunt for a defect in `helicoid` — the corpus is the instrument,
and `0006` puts the instrument first precisely so that this class of error is not discovered through
the thing being measured. Generating it right costs one afternoon in the generator; generating it
wrong costs an unknown amount of time spent not believing a correct implementation.

For §3 the alternative is to leave the wording alone and let the implementer pick the law. One of
the two choices makes a NORMATIVE test assert a falsehood, and the other needs a fixture inequality
nothing in the spec states. §3's instruction — "do not fix it" — is exactly right about the
mathematics and silent about the one thing an implementer has to decide, which is the kind of gap
`CLAUDE.md` says to stop and ask about rather than guess.

## Consequences

- `PHASE4.md` §3 and §4 are edited; both stay NORMATIVE and both become executable.
- The generator gains a geodesic module whose reference does not use `mp.logm`, which is the first
  corpus id in the repository to depart from "`mp.expm`/`mp.logm` of the matrix". The same departure
  is owed by `so3_log`, `se2_log`, `s2_local` and `sim3_log`; this record does not take it for them,
  and the `docs/maths/index.md` items stay open.
- The near-π rows of both ids are documented as conditioning measurements. A domination comparison
  on them is a comparison of two sign conventions as much as of two kernels, which `PHASE4.md` §5.2
  should read with that in mind.
- `PHASE1.md` §4.3 gains two rows and the `geo:*` family, in the same step as the corpus files
  (the table is machine-equated to `coverage::required`); the corpus grows by two files, from
  34 MB of the 50 MB budget.

## Implementation plan

1. This record and `PHASE4.md` §3 and §4's edits — verified by `just lint`. **Not** the
   `PHASE1.md` §4.3 rows: `envelope/coverage.rs`'s
   `the_required_ids_are_the_ids_of_the_definitions_table` equates §4.3's ids with
   `coverage::required` exactly, and both geodesic ids are in `coverage::OWED` until their corpus
   files exist, so a §4.3 row written early fails that test. The row, the `OWED`-to-`required`
   move and the corpus file are one step, below.
2. `gen/strata.py`'s `geo:*` family, `gen/geodesic.py`, `gen/check_geodesic.py`, the
   `gen/registry.py` entries, and the generator's own prefix test extended — verified by
   `just corpus-test`.
3. `just corpus`, the committed files, the two `PHASE1.md` §4.3 definition rows and the
   `coverage::OWED`-to-`coverage::required` move, with the §4 cross-check green at
   $\theta \le 3.0$ and the identity checks green everywhere — verified by `just corpus-check`
   byte-comparing a regeneration and by `the_required_ids_are_the_ids_of_the_definitions_table`.
4. §3's three invariance tests with their measured bounds, the right-invariance fixture carrying
   GE.5(b)'s inequality in its comment — verified by `just test`.

## Open questions

None.

## Further work

1. `so3_log`'s own `mp.logm` onset ($\theta \approx 3.03$) is the same defect in a committed corpus
   id, and `so3_log`'s `theta:pi-1e-k` strata are inside it. Whether those references are wrong
   today is a question this record does not answer and the `docs/maths/index.md` item holds.
2. A corpus id for the geodesic **Jacobians**; §4 has none, and
   [`0043`](./0043-the-geodesic-jacobian-ships-the-cancellation-free-form.md) notes it too.
3. `geo:consecutive` overlaps `tf_tree`'s own `interp-accuracy` fixture. Whether the two should be
   the same inputs — so that a `tf_tree` regression and a `helicoid` stratum name the same cell —
   is worth asking when Wave 3 lands.
