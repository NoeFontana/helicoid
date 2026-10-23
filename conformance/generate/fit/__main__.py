import argparse
from pathlib import Path

from . import atan2


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="fit", description="candidate kernels' polynomials")
    sub = parser.add_subparsers(dest="command", required=True)
    one = sub.add_parser("atan2", help="the portable atan2 candidate (docs/decisions/0105, draft)")
    one.add_argument("--f64", required=True, metavar="CENTRES:TERMS", help="e.g. 8:5")
    one.add_argument("--f32", required=True, metavar="TERMS", type=int, help="e.g. 8")
    one.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    centres, terms = (int(v) for v in args.f64.split(":"))
    args.out.write_text(atan2.write(centres, terms, args.f32))
    print(args.out.read_text(), end="")
    return 0


raise SystemExit(main())
