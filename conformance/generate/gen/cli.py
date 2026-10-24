import argparse
import os
from collections import Counter
from pathlib import Path

from . import corpus, manifest, memo
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
    everything.add_argument(
        "--memo", type=Path, help="replay unchanged tasks from this directory (docs/decisions/0108)"
    )
    everything.add_argument(
        "--fresh",
        action="store_true",
        help="with --memo, run every task and write its entry, replaying none",
    )
    everything.add_argument(
        "--audit",
        type=float,
        default=0.02,
        help="with --memo, the share of hits run again and compared (default 0.02)",
    )
    verify = sub.add_parser(
        "verify", help="a corpus's manifest names this generator and every file's hash"
    )
    verify.add_argument("--corpus", type=Path, required=True)
    key = sub.add_parser(
        "key", help="what a full check of this corpus by this generator proves, as one hash"
    )
    key.add_argument("--corpus", type=Path, required=True)
    deep = sub.add_parser(
        "deep", help="the deep tier, with HELICOID_CORPUS_DEEP=1 (docs/decisions/0088)"
    )
    deep.add_argument("--out", type=Path, required=True)
    deep.add_argument("--committed", type=Path, required=True, help="the committed corpus")
    deep.add_argument("--jobs", type=int, help="worker processes (default: every core)")
    args = parser.parse_args(argv)
    if args.command == "list":
        setup()
        for name, spec in FUNCTIONS.items():
            n = sum(len(spec.inputs(s)) for s in spec.strata)
            print(f"{name}\t{n} records\t{len(spec.strata)} strata")
    elif args.command == "deep":
        for name, built in corpus.write_deep(args.out, args.committed, args.jobs).items():
            print(f"{name}: {built.records} records, {built.checked} {CHECKED[built.kind]}")
    elif args.command == "verify":
        problems = manifest.verify(args.corpus)
        for problem in problems:
            print(f"corpus-verify: {problem}")
        if problems:
            return 1
        print(f"corpus-verify: {args.corpus} is what its manifest says, for this generator")
    elif args.command == "key":
        print(manifest.check_key(args.corpus))
    else:
        cache = None
        if args.memo is not None:
            seed = os.environ.get("GITHUB_RUN_ID") or os.urandom(8).hex()
            cache = memo.Memo(str(args.memo.resolve()), args.audit, seed, replay=not args.fresh)
        stats: Counter = Counter()
        for name, built in corpus.write(args.out, args.jobs, cache, stats).items():
            print(f"{name}: {built.records} records, {built.checked} {CHECKED[built.kind]}")
        if cache is not None:
            print(
                f"memo: {stats['hit']} tasks replayed, {stats['audited']} replayed and run again, "
                f"{stats['fresh']} run"
            )
    return 0
