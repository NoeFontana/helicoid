"""Corpus assembly: one JSONL file per function id, records in canonical stratum order."""

import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from mpmath import mp, mpf

from . import fmt, manifest, series, switchref
from .fmt import Mat
from .manifest import Built
from .precision import RECHECK_DIGITS, RECHECK_DPS, is_binary32, setup
from .registry import FUNCTIONS, FunctionSpec
from .strata import Stratum


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


def recheck(spec: FunctionSpec, inputs: dict, out: dict) -> None:
    with mp.workdps(RECHECK_DPS):
        hi = spec.evaluate(inputs)
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


def recheck_stratum(spec: FunctionSpec, inputs: list[dict], outs: list[dict]) -> int:
    """Recheck the picked records of one stratum; how many is what the manifest reports."""
    picks = sorted(recheck_indices(len(inputs)))
    for i in picks:
        recheck(spec, inputs[i], outs[i])
    return len(picks)


def _encode(values: dict, one) -> dict:
    """A vector is an array; a `Mat` is its column-major array and a sibling `"shape"` key, so a
    dict holds at most one."""
    out: dict = {}
    for k, v in values.items():
        if isinstance(v, Mat):
            if "shape" in out:
                raise ValueError(f"two matrices in one record part: {sorted(values)}")
            out[k], out["shape"] = [one(x) for x in v.data], list(v.shape)
        else:
            out[k] = [one(x) for x in v] if isinstance(v, (list, tuple)) else one(v)
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
    if stratum.f32:
        require_binary32(spec, stratum, inputs)
    outs = [spec.evaluate(inp) for inp in inputs]
    rechecked = recheck_stratum(spec, inputs, outs)
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


def write(out_dir: Path, jobs: int | None = None) -> dict[str, Built]:
    """Write every file and the manifest; a `*.jsonl` no function id owns any more is removed."""
    out_dir.mkdir(parents=True, exist_ok=True)
    files = {f"{name}.jsonl": Built(*built) for name, built in build_all(jobs).items()}
    files["coeff_series.jsonl"] = Built(*series.build(), kind="series")
    files["coeff_switch_ref.jsonl"] = Built(*switchref.build(), kind="switch-ref")
    for stale in out_dir.glob("*.jsonl"):
        if stale.name not in files:
            stale.unlink()
    for name, built in files.items():
        (out_dir / name).write_bytes(built.data)
    (out_dir / "MANIFEST.json").write_bytes(manifest.render(files))
    return files
