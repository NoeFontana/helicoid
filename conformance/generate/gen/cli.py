import argparse
from pathlib import Path

from . import corpus
from .manifest import CHECKED
from .precision import setup
from .registry import FUNCTIONS


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="gen", description="helicoid conformance corpus generator"
    )
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="function ids and their record counts")
    everything = sub.add_parser("all", help="generate every file")
    everything.add_argument("--out", type=Path, required=True)
    everything.add_argument("--jobs", type=int, help="worker processes (default: every core)")
    args = parser.parse_args(argv)
    if args.command == "list":
        setup()
        for name, spec in FUNCTIONS.items():
            n = sum(len(spec.inputs(s)) for s in spec.strata)
            print(f"{name}\t{n} records\t{len(spec.strata)} strata")
    else:
        for name, built in corpus.write(args.out, args.jobs).items():
            print(f"{name}: {built.records} records, {built.checked} {CHECKED[built.kind]}")
    return 0
