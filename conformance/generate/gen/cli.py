import argparse
from pathlib import Path

from . import corpus
from .precision import RECHECK_DPS
from .registry import FUNCTIONS


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="gen", description="helicoid conformance corpus generator"
    )
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="function ids and their record counts")
    sub.add_parser("all", help="generate every file").add_argument(
        "--out", type=Path, required=True
    )
    args = parser.parse_args(argv)
    if args.command == "list":
        for name, spec in FUNCTIONS.items():
            print(f"{name}\t{sum(s.count for s in spec.strata)} records\t{len(spec.strata)} strata")
    else:
        for name, (_, records, rechecked) in corpus.write(args.out).items():
            print(f"{name}: {records} records, {rechecked} rechecked at {RECHECK_DPS} digits")
    return 0
