"""MANIFEST.json: what the corpus was generated with, and a hash of every file in it."""

import hashlib
import json
import sys
from pathlib import Path
from typing import NamedTuple

import mpmath

from .precision import DPS, RECHECK_DIGITS, RECHECK_DPS
from .rng import SEED

ROOT = Path(__file__).resolve().parent.parent
# What `Built.checked` counts, per file kind: "corpus" files hold (id, stratum, in, out) records,
# recomputed at 150 digits; "series" (`coeff_series`) holds series, each equal to exact algebra.
CHECKED = {"corpus": "rechecked", "series": "verified"}


class Built(NamedTuple):
    data: bytes
    records: int
    checked: int
    kind: str = "corpus"


def generator_identity(root: Path = ROOT) -> str:
    """SHA-256 over the sorted (path, bytes) of everything that determines the corpus.

    A content hash, not a git revision: a commit hash would change with every commit and break
    `corpus-check`, and a commit cannot contain its own hash.
    """
    files = [root / "pyproject.toml", root / "uv.lock", root / ".python-version"]
    files += (root / "gen").rglob("*.py")
    h = hashlib.sha256()
    for path in sorted(files, key=lambda p: p.relative_to(root).as_posix()):
        data = path.read_bytes()
        h.update(f"{path.relative_to(root).as_posix()}\0{len(data)}\0".encode() + data + b"\n")
    return "sha256:" + h.hexdigest()


def render(files: dict[str, Built], root: Path = ROOT) -> bytes:
    """`files` maps a file name to what was built; each entry says its `kind` and `CHECKED[kind]`."""
    manifest = {
        "dps": DPS,
        "files": {
            name: {
                "kind": b.kind,
                CHECKED[b.kind]: b.checked,
                "records": b.records,
                "sha256": hashlib.sha256(b.data).hexdigest(),
            }
            for name, b in sorted(files.items())
        },
        "generator": generator_identity(root),
        "mpmath": mpmath.__version__,
        "python": f"{sys.version_info.major}.{sys.version_info.minor}",  # never the micro version
        "recheck": {"digits": RECHECK_DIGITS, "dps": RECHECK_DPS},
        "seed": f"0x{SEED:016x}",
    }
    return (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode()
