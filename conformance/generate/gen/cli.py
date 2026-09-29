import argparse
from pathlib import Path

from . import corpus
from .manifest import CHECKED
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
        for name, built in corpus.write(args.out).items():
            print(f"{name}: {built.records} records, {built.checked} {CHECKED[built.kind]}")
    return 0
