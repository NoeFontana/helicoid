"""Memoised corpus tasks (docs/decisions/0108).

A task is one (id, stratum) of the registry, or one of the files `corpus.write` builds beside it.
Its output is a function of four things, and the memo replays it only while all four are unchanged:

- the code it ran: every file of `gen` in which a function started during the task
  (`sys.monitoring`, one event per code object per task);
- the module-level values that code read: the static import closure of those files within `gen`,
  without `registry` or `cli`, whose only effect on a task is the spec it names (unless their own
  code ran);
- the spec or builder that names it: `describe` of it, by qualified name, code, closure values and
  fields;
- the interpreter: Python's full version, mpmath's, the platform and HELICOID_CORPUS_DEEP.

A hit is re-run with probability `audit`, and a byte that differs fails the run: the memo is an
optimisation of `corpus-check`, never a substitute for the fresh run that `main` and the nightly
job make.
"""

import ast
import dataclasses
import functools
import hashlib
import json
import marshal
import os
import platform
import sys
import types
from collections.abc import Callable
from pathlib import Path
from typing import Any

import mpmath

GEN = Path(__file__).resolve().parent
# Their effect on a task is the spec `describe` reads: an import of them is not a dependency.
NOT_IMPORTED = frozenset({"registry.py", "cli.py", "__main__.py"})
# Bumped when an entry's layout or meaning changes.
FORMAT = 1


class MemoError(Exception):
    pass


@dataclasses.dataclass(frozen=True)
class Memo:
    """Where the entries live, and the share of hits re-run (`audit`, keyed by `seed`). With
    `replay` false every task runs and its entry is written: a fresh run that fills the memo."""

    root: str
    audit: float = 0.0
    seed: str = ""
    replay: bool = True


def environment() -> str:
    deep = os.environ.get("HELICOID_CORPUS_DEEP", "")
    return f"{FORMAT}|{sys.version}|{mpmath.__version__}|{platform.machine()}|{deep}"


def _canonical(code: types.CodeType) -> types.CodeType:
    """`code` without its file name, nested code included, so the same source marshals to the
    same bytes in every checkout."""
    consts = tuple(_canonical(c) if isinstance(c, types.CodeType) else c for c in code.co_consts)
    return code.replace(co_filename="", co_consts=consts)


def describe(x: Any, seen: frozenset[int] = frozenset()) -> str:
    """A text that changes whenever `x` could make a task compute something else. An object it
    cannot describe raises TypeError, and that task is not memoised."""
    if x is None or isinstance(x, (bool, int, float, str, bytes)):
        return repr(x)
    if isinstance(x, mpmath.mpf):
        return f"mpf{x._mpf_}"  # exact, whatever the working precision
    if id(x) in seen:
        return "<cycle>"
    seen = seen | {id(x)}
    if isinstance(x, (tuple, list)):
        return "[" + ",".join(describe(v, seen) for v in x) + "]"
    if isinstance(x, (set, frozenset)):
        return "{" + ",".join(sorted(describe(v, seen) for v in x)) + "}"
    if isinstance(x, dict):
        items = sorted((describe(k, seen), describe(v, seen)) for k, v in x.items())
        return "{" + ",".join(f"{k}:{v}" for k, v in items) + "}"
    if isinstance(x, functools.partial):
        return f"partial({describe(x.func, seen)},{describe(x.args, seen)},{describe(x.keywords, seen)})"
    if isinstance(x, types.MethodType):
        return f"method({describe(x.__func__, seen)},{describe(x.__self__, seen)})"
    if isinstance(x, types.FunctionType):
        cells = []
        for cell in x.__closure__ or ():
            try:
                cells.append(describe(cell.cell_contents, seen))
            except ValueError:
                cells.append("<empty>")
        digest = hashlib.sha256(marshal.dumps(_canonical(x.__code__))).hexdigest()
        return (
            f"fn({x.__module__}.{x.__qualname__},{digest},[{','.join(cells)}],"
            f"{describe(x.__defaults__, seen)},{describe(x.__kwdefaults__, seen)})"
        )
    if isinstance(x, type):
        return f"type({x.__module__}.{x.__qualname__})"
    if isinstance(x, (types.BuiltinFunctionType, types.ModuleType)):
        return f"{type(x).__name__}({getattr(x, '__module__', '')}.{getattr(x, '__qualname__', x.__name__)})"
    if dataclasses.is_dataclass(x):
        fields = ",".join(
            f"{f.name}={describe(getattr(x, f.name), seen)}" for f in dataclasses.fields(x)
        )
        return f"{type(x).__qualname__}({fields})"
    raise TypeError(f"memo: cannot describe a {type(x).__qualname__}")


# --- what a task ran ----------------------------------------------------------------------------

_ran: set[str] = set()
_tool: int | None = None


def _started(code: types.CodeType, _offset: int) -> Any:
    _ran.add(code.co_filename)
    return sys.monitoring.DISABLE


def _watch() -> None:
    """Start recording which files' functions run; once per process, then re-armed per task."""
    global _tool
    if _tool is None:
        mon = sys.monitoring
        _tool = next(i for i in range(6) if mon.get_tool(i) is None)
        mon.use_tool_id(_tool, "helicoid-corpus-memo")
        mon.register_callback(_tool, mon.events.PY_START, _started)
        mon.set_events(_tool, mon.events.PY_START)
    _ran.clear()
    sys.monitoring.restart_events()


@functools.cache
def _imports(path: Path) -> frozenset[Path]:
    """The modules of `gen` that `path` imports, statically."""
    names: set[str] = set()
    for node in ast.walk(ast.parse(path.read_text(), str(path))):
        if isinstance(node, ast.ImportFrom):
            if node.level == 1 and node.module:
                names.add(node.module.split(".")[0])
            elif node.level == 1:
                names.update(a.name for a in node.names)
            elif node.level == 0 and node.module and node.module.split(".")[0] == "gen":
                parts = node.module.split(".")
                names.update([parts[1]] if len(parts) > 1 else [a.name for a in node.names])
        elif isinstance(node, ast.Import):
            for a in node.names:
                parts = a.name.split(".")
                if parts[0] == "gen" and len(parts) > 1:
                    names.add(parts[1])
    return frozenset(GEN / f"{n}.py" for n in names if (GEN / f"{n}.py").is_file())


def dependencies(ran: set[str]) -> list[Path]:
    """The files of `gen` a task depends on: those it ran, `__init__`, and their import closure
    without `NOT_IMPORTED`."""
    todo = [Path(f) for f in ran if Path(f).parent == GEN] + [GEN / "__init__.py"]
    found: set[Path] = set()
    while todo:
        path = todo.pop()
        if path in found:
            continue
        found.add(path)
        todo.extend(p for p in _imports(path) if p.name not in NOT_IMPORTED)
    return sorted(found)


@functools.cache
def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# --- the memo -----------------------------------------------------------------------------------


def _entry(memo: Memo, task: tuple[str, str]) -> Path:
    name = hashlib.sha256("\0".join(task).encode()).hexdigest()[:32]
    return Path(memo.root) / f"{name}.json"


def _audited(memo: Memo, key: str) -> bool:
    if memo.audit <= 0:
        return False
    draw = int(hashlib.sha256(f"{memo.seed}\0{key}".encode()).hexdigest()[:8], 16) / 2**32
    return draw < memo.audit


def run(memo: Memo | None, task: tuple[str, str], spec: Any, compute: Callable[[], Any]) -> tuple[Any, str]:
    """`compute()`, or its memoised value, and which: `fresh`, `hit` or `audited`. `compute`
    returns a JSON value; `spec` is what names the task (see `describe`)."""
    if memo is None:
        return compute(), "fresh"
    try:
        key = hashlib.sha256(
            f"{environment()}\0{task}\0{describe(spec)}".encode()
        ).hexdigest()
    except TypeError:
        return compute(), "fresh"
    path = _entry(memo, task)
    try:
        entry = json.loads(path.read_text())
    except (OSError, ValueError):
        entry = None
    if memo.replay and entry is not None and entry["key"] == key:
        deps = entry["deps"]
        if all((GEN / rel).is_file() and _hash(GEN / rel) == h for rel, h in deps.items()):
            if not _audited(memo, key):
                return entry["value"], "hit"
            value = _fresh(memo, task, key, path, compute)
            if json.dumps(value, sort_keys=True) != json.dumps(entry["value"], sort_keys=True):
                raise MemoError(
                    f"memo: {task} replayed a value its fresh run does not give; "
                    "rerun with --fresh and report it (docs/decisions/0108)"
                )
            return value, "audited"
    return _fresh(memo, task, key, path, compute), "fresh"


def _fresh(memo: Memo, task: tuple[str, str], key: str, path: Path, compute: Callable[[], Any]) -> Any:
    _watch()
    value = compute()
    ran = set(_ran)
    deps = {p.relative_to(GEN).as_posix(): _hash(p) for p in dependencies(ran)}
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(f".{os.getpid()}.tmp")
    tmp.write_text(json.dumps({"task": list(task), "key": key, "deps": deps, "value": value}))
    tmp.replace(path)
    return value
