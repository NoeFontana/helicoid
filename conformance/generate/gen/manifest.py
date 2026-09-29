"""MANIFEST.json: what the corpus was generated with, and a hash of every file in it."""

import hashlib
import json
import sys
from pathlib import Path

import mpmath

from .precision import DPS, RECHECK_DIGITS, RECHECK_DPS
from .rng import SEED

ROOT = Path(__file__).resolve().parent.parent


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


def render(files: dict[str, tuple[bytes, int, int]], root: Path = ROOT) -> bytes:
    """`files` maps a file name to (bytes, record count, rechecked count)."""
    manifest = {
        "dps": DPS,
        "files": {
            name: {"records": n, "rechecked": r, "sha256": hashlib.sha256(data).hexdigest()}
            for name, (data, n, r) in sorted(files.items())
        },
        "generator": generator_identity(root),
        "mpmath": mpmath.__version__,
        "python": f"{sys.version_info.major}.{sys.version_info.minor}",  # never the micro version
        "recheck": {"digits": RECHECK_DIGITS, "dps": RECHECK_DPS},
        "seed": f"0x{SEED:016x}",
    }
    return (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode()
