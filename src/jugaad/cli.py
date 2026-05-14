"""Top-level command line interface for jugaad."""

from __future__ import annotations

import argparse
import signal
import sys
from collections.abc import Sequence

from jugaad import __version__
from jugaad.commands import find_retained_introns, gff3_to_introns


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="jugaad",
        description="Small bioinformatics command line utilities.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    subparsers = parser.add_subparsers(dest="command", metavar="<command>")
    find_retained_introns.register_parser(subparsers)
    gff3_to_introns.register_parser(subparsers)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    if hasattr(signal, "SIGPIPE"):
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    parser = build_parser()
    args_list = list(sys.argv[1:] if argv is None else argv)
    if not args_list:
        parser.print_help(sys.stderr)
        return 1

    args = parser.parse_args(args_list)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
