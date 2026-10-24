"""The generator's unit tests, one test class per task across every core (docs/decisions/0108).

`python -m tests.parallel [--jobs N]` from conformance/generate. A class is the unit, so a
`setUpClass` runs once as it does serially; a class's tests run in their written order. The exit
status and the failures are `unittest`'s.
"""

import argparse
import io
import os
import sys
import time
import unittest
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent


def classes() -> list[str]:
    """Every test class's dotted name, the largest modules first so the pool starts them first."""
    suite = unittest.defaultTestLoader.discover(str(HERE))
    names: dict[str, None] = {}

    def walk(s):
        for t in s:
            if isinstance(t, unittest.TestSuite):
                walk(t)
            elif isinstance(t, unittest.loader._FailedTest):
                raise SystemExit(f"tests.parallel: {t.id()} does not import:\n{t._exception}")
            else:
                names[f"{type(t).__module__}.{type(t).__qualname__}"] = None

    walk(suite)
    size = lambda name: (HERE / f"{name.split('.')[0]}.py").stat().st_size
    return sorted(names, key=lambda n: (-size(n), n))


def run(name: str) -> tuple[str, bool, int, int, str, float]:
    """One class: (name, passed, tests run, skipped, unittest's report, seconds)."""
    stream = io.StringIO()
    start = time.perf_counter()
    suite = unittest.defaultTestLoader.loadTestsFromName(name)
    result = unittest.TextTestRunner(stream=stream, verbosity=0).run(suite)
    took = time.perf_counter() - start
    return name, result.wasSuccessful(), result.testsRun, len(result.skipped), stream.getvalue(), took


def main() -> int:
    parser = argparse.ArgumentParser(prog="tests.parallel")
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    args = parser.parse_args()
    # As `unittest discover -s tests` from conformance/generate: `gen` from here, tests by module.
    sys.path[:0] = [str(HERE.parent), str(HERE)]
    start = time.perf_counter()
    names = classes()
    with ProcessPoolExecutor(args.jobs) as pool:
        results = list(pool.map(run, names))
    failed = [r for r in results if not r[1]]
    for name, _, _, _, report, _ in failed:
        print(f"=== {name}\n{report}")
    tests, skipped = sum(r[2] for r in results), sum(r[3] for r in results)
    slowest = max(results, key=lambda r: r[5])
    print(
        f"Ran {tests} tests ({skipped} skipped) in {len(names)} classes in "
        f"{time.perf_counter() - start:.1f}s on {args.jobs} workers; slowest {slowest[0]} "
        f"{slowest[5]:.1f}s"
    )
    print("FAILED" if failed else "OK")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
