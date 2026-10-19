"""Corpus assembly: one JSONL file per function id, records in canonical stratum order."""

import json

import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from mpmath import mp, mpf

from . import fmt, manifest, series, sim3coeff, switchref
from .fmt import Mat
from .manifest import Built
from .precision import RECHECK_DIGITS, RECHECK_DPS, is_binary32, setup
from .registry import FUNCTIONS, FunctionSpec
from .strata import DEEP, Stratum


class RecheckError(Exception):
    pass


class BinaryError(Exception):
    pass


def _flat(value) -> list:
    if isinstance(value, Mat):
        return list(value.data)
    return list(value) if isinstance(value, (list, tuple)) else [value]


def recheck_indices(n: int) -> set[int]:
    """Evenly spaced, at least 1% of a stratum and at least one record of every stratum."""
    m = max(1, -(-n // 100))
    return {(2 * j + 1) * n // (2 * m) for j in range(m)}


def evaluate(spec: FunctionSpec, inputs: dict, f32: bool) -> dict:
    """The reference: of the inputs alone, or for a conversion (0087) at the precision too."""
    return spec.evaluate(inputs, f32) if spec.converts else spec.evaluate(inputs)


def recheck(spec: FunctionSpec, inputs: dict, out: dict, f32: bool = False) -> None:
    with mp.workdps(RECHECK_DPS):
        hi = evaluate(spec, inputs, f32)
    if out.keys() != hi.keys():
        raise RecheckError(f"{spec.name} {inputs}: outputs {sorted(out)} vs {sorted(hi)}")
    tol = mpf(10) ** -RECHECK_DIGITS
    for key, value in out.items():
        lo_values, hi_values = _flat(value), _flat(hi[key])
        if len(lo_values) != len(hi_values):
            raise RecheckError(f"{spec.name} {key} {inputs}: {len(lo_values)} vs {len(hi_values)}")
        for a, b in zip(lo_values, hi_values, strict=True):
            if abs(a - b) > tol * max(abs(a), abs(b)):
                raise RecheckError(f"{spec.name} {key} {inputs}: dps disagreement, {a} vs {b}")


def recheck_stratum(
    spec: FunctionSpec, inputs: list[dict], outs: list[dict], f32: bool = False
) -> int:
    """Recheck the picked records of one stratum; how many is what the manifest reports."""
    picks = sorted(recheck_indices(len(inputs)))
    for i in picks:
        recheck(spec, inputs[i], outs[i], f32)
    return len(picks)


def _encode(values: dict, one) -> dict:
    """A vector is an array; a `Mat` is its column-major array beside a `"shape"` key: `[r, c]`
    where the part holds one `Mat`, so no record of one moves, and `{field: [r, c], ...}` where it
    holds several (docs/decisions/0101 decision 8)."""
    if "shape" in values:
        raise ValueError(f"a field named `shape` in a record part: {sorted(values)}")
    out: dict = {}
    shapes: dict = {}
    for k, v in values.items():
        if isinstance(v, Mat):
            out[k], shapes[k] = [one(x) for x in v.data], list(v.shape)
        else:
            out[k] = [one(x) for x in v] if isinstance(v, (list, tuple)) else one(v)
    if shapes:
        out["shape"] = next(iter(shapes.values())) if len(shapes) == 1 else shapes
    return out


def require_binary32(spec: FunctionSpec, stratum: Stratum, inputs: list[dict]) -> None:
    """Every input of an `@f32` stratum is exactly a binary32, so a binary32 subject receives it
    by a lossless cast and the reference is the function at what it evaluates (0016)."""
    for inp in inputs:
        for key, value in inp.items():
            for x in _flat(value):
                if not is_binary32(x):
                    raise BinaryError(f"{spec.name} {stratum.name}: {key} = {x.hex()}")


def _build_stratum(spec: FunctionSpec, stratum: Stratum) -> tuple[list[dict], int]:
    """The records of one stratum (without their ids) and how many were rechecked."""
    setup()
    inputs = spec.inputs(stratum)
    if stratum.f32 and not spec.converts:
        require_binary32(spec, stratum, inputs)
    outs = [evaluate(spec, inp, stratum.f32) for inp in inputs]
    rechecked = recheck_stratum(spec, inputs, outs, stratum.f32)
    if spec.check is not None:
        for inp, out in zip(inputs, outs, strict=True):
            spec.check(inp, out)
    records = [
        {"stratum": stratum.name, "in": _encode(inp, fmt.hex_float), "out": _encode(out, fmt.dec30)}
        for inp, out in zip(inputs, outs, strict=True)
    ]
    return records, rechecked


def _assemble(parts: list[tuple[list[dict], int]]) -> tuple[bytes, int, int]:
    lines = []
    for records, _ in parts:
        lines += [fmt.dumps({"id": len(lines) + i, **r}) + "\n" for i, r in enumerate(records)]
    return "".join(lines).encode(), len(lines), sum(n for _, n in parts)


def build(spec: FunctionSpec, only: set[str] | None = None) -> tuple[bytes, int, int]:
    """(file bytes, record count, rechecked count). `only` restricts strata (tests)."""
    return _assemble(
        [_build_stratum(spec, s) for s in spec.strata if only is None or s.name in only]
    )


def _registered(name: str, stratum: str) -> tuple[list[dict], int]:
    """A worker's task: names only, since a spec holds lambdas and does not pickle."""
    spec = FUNCTIONS[name]
    return _build_stratum(spec, next(s for s in spec.strata if s.name == stratum))


def build_all(jobs: int | None = None) -> dict[str, tuple[bytes, int, int]]:
    """Every registered id. One task per (id, stratum), largest first, joined in catalogue order:
    the bytes are those of a serial run whatever `jobs` (default: every core)."""
    tasks = [(name, s.name, s.count) for name, spec in FUNCTIONS.items() for s in spec.strata]
    jobs = jobs or os.cpu_count() or 1
    if jobs == 1:
        results = [_registered(name, stratum) for name, stratum, _ in tasks]
    else:
        with ProcessPoolExecutor(jobs) as pool:
            futures = {
                i: pool.submit(_registered, name, stratum)
                for i, (name, stratum, count) in sorted(enumerate(tasks), key=lambda t: -t[1][2])
            }
            results = [futures[i].result() for i in range(len(tasks))]
    parts: dict[str, list] = {name: [] for name in FUNCTIONS}
    for (name, _, _), result in zip(tasks, results, strict=True):
        parts[name].append(result)
    return {name: _assemble(p) for name, p in parts.items()}


class PrefixError(Exception):
    pass


def _by_stratum(data: bytes) -> dict[str, list[str]]:
    """A file's records per stratum, in order, each without its `id`: what nesting compares."""
    out: dict[str, list[str]] = {}
    for line in data.decode().splitlines():
        record = json.loads(line)
        record.pop("id")
        out.setdefault(record["stratum"], []).append(fmt.dumps(record))
    return out


def check_prefix(name: str, committed: bytes, deep: bytes) -> None:
    """Every committed stratum of `name` is the first records of its deep namesake, byte for byte
    but the id, and the deep file has no stratum the committed one lacks (docs/decisions/0088)."""
    wide, narrow = _by_stratum(deep), _by_stratum(committed)
    if list(wide) != list(narrow):
        raise PrefixError(f"{name}: strata {list(narrow)} committed, {list(wide)} deep")
    for stratum, records in narrow.items():
        if wide[stratum][: len(records)] != records:
            raise PrefixError(f"{name} {stratum}: the committed records are not the deep prefix")


def write_deep(out_dir: Path, committed: Path, jobs: int | None = None) -> dict[str, Built]:
    """The deep tier (docs/decisions/0088): every id whose records grow under
    HELICOID_CORPUS_DEEP=1, each checked against the committed corpus with [`check_prefix`], and
    the manifest of those files alone. The series and switch references do not grow."""
    if not DEEP:
        raise SystemExit("gen deep draws the deep tier: run it with HELICOID_CORPUS_DEEP=1")
    listed = json.loads((committed / "MANIFEST.json").read_text())["files"]
    files: dict[str, Built] = {}
    for name, built in build_all(jobs).items():
        file = f"{name}.jsonl"
        if built[1] == listed[file]["records"]:
            continue
        check_prefix(name, (committed / file).read_bytes(), built[0])
        files[file] = Built(*built)
    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in out_dir.glob("*.jsonl"):
        if stale.name not in files:
            stale.unlink()
    for name, built in files.items():
        (out_dir / name).write_bytes(built.data)
    (out_dir / "MANIFEST.json").write_bytes(manifest.render(files))
    return files


def write(out_dir: Path, jobs: int | None = None) -> dict[str, Built]:
    """Write every file and the manifest; a `*.jsonl` no function id owns any more is removed.
    Refused under HELICOID_CORPUS_DEEP=1, whose sample counts would rewrite the committed corpus at
    the deep tier's size (docs/decisions/0088)."""
    if DEEP:
        raise SystemExit(
            "gen all writes the committed corpus: unset HELICOID_CORPUS_DEEP, or run `gen deep`"
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    files = {f"{name}.jsonl": Built(*built) for name, built in build_all(jobs).items()}
    files["coeff_series.jsonl"] = Built(*series.build(), kind="series")
    files["coeff_switch_ref.jsonl"] = Built(*switchref.build(), kind="switch-ref")
    files["sim3_series.jsonl"] = Built(*sim3coeff.build_series(), kind="series")
    files["sim3_ref.jsonl"] = Built(*sim3coeff.build_ref(), kind="switch-ref")
    for stale in out_dir.glob("*.jsonl"):
        if stale.name not in files:
            stale.unlink()
    for name, built in files.items():
        (out_dir / name).write_bytes(built.data)
    (out_dir / "MANIFEST.json").write_bytes(manifest.render(files))
    return files
