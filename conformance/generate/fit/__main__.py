import argparse
from fractions import Fraction
from pathlib import Path

from . import atan2


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="fit", description="candidate kernels' polynomials")
    sub = parser.add_subparsers(dest="command", required=True)
    one = sub.add_parser("atan2", help="the portable atan2 candidate (docs/decisions/0105, draft)")
    one.add_argument(
        "--fit",
        action="append",
        required=True,
        metavar="PRECISION:TERMS:UMAX",
        help="one per precision, e.g. f64:13:1/2",
    )
    one.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    specs = []
    for spec in args.fit:
        precision, terms, umax = spec.split(":")
        specs.append((precision, int(terms), Fraction(umax)))
    args.out.write_text(atan2.write(specs))
    print(args.out.read_text(), end="")
    return 0


raise SystemExit(main())
