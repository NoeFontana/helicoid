"""Corpus assembly: one JSONL file per function id, records in canonical stratum order."""

from pathlib import Path

from mpmath import mp, mpf

from . import fmt, manifest, series
from .manifest import Built
from .precision import RECHECK_DIGITS, RECHECK_DPS, setup
from .registry import FUNCTIONS, FunctionSpec


class RecheckError(Exception):
    pass


def _flat(value) -> list:
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
    return {
        k: [one(x) for x in v] if isinstance(v, (list, tuple)) else one(v)
        for k, v in values.items()
    }


def build(spec: FunctionSpec, only: set[str] | None = None) -> tuple[bytes, int, int]:
    """(file bytes, record count, rechecked count). `only` restricts strata (tests)."""
    setup()
    lines, rechecked = [], 0
    for stratum in spec.strata:
        if only is not None and stratum.name not in only:
            continue
        inputs = spec.inputs(stratum)
        outs = [spec.evaluate(inp) for inp in inputs]
        rechecked += recheck_stratum(spec, inputs, outs)
        if spec.check is not None:
            for inp, out in zip(inputs, outs, strict=True):
                spec.check(inp, out)
        for inp, out in zip(inputs, outs, strict=True):
            record = {
                "id": len(lines),
                "stratum": stratum.name,
                "in": _encode(inp, fmt.hex_float),
                "out": _encode(out, fmt.dec30),
            }
            lines.append(fmt.dumps(record) + "\n")
    return "".join(lines).encode(), len(lines), rechecked


def write(out_dir: Path) -> dict[str, Built]:
    """Write every file and the manifest; a `*.jsonl` no function id owns any more is removed."""
    out_dir.mkdir(parents=True, exist_ok=True)
    files = {f"{name}.jsonl": Built(*build(spec)) for name, spec in FUNCTIONS.items()}
    files["coeff_series.jsonl"] = Built(*series.build(), kind="series")
    for stale in out_dir.glob("*.jsonl"):
        if stale.name not in files:
            stale.unlink()
    for name, built in files.items():
        (out_dir / name).write_bytes(built.data)
    (out_dir / "MANIFEST.json").write_bytes(manifest.render(files))
    return files
