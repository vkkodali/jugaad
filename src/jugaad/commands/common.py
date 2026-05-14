"""Shared helpers for command modules."""

from __future__ import annotations

import argparse
import gzip
import sys
from contextlib import ExitStack
from datetime import datetime
from typing import TextIO


def currtime() -> str:
    return datetime.now().isoformat(timespec="seconds", sep=" ")


def add_gxf_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("-i", "--infile", help="input file; default STDIN")
    parser.add_argument("-o", "--outfile", help="output file; default STDOUT")
    parser.add_argument(
        "-z",
        "--gzip",
        help="input is gzip compressed",
        action="store_true",
    )
    parser.add_argument(
        "-f",
        "--format",
        help="input file type",
        choices=["gtf", "gff3"],
        required=True,
    )
    parser.add_argument(
        "--tx_attrib",
        help="attribute to use for transcript id; default is `transcript_id`",
        default="transcript_id",
    )
    parser.add_argument(
        "--gene_attrib",
        help="attribute to use for gene id; default is `GeneID`",
        default="GeneID",
    )


def open_input(args: argparse.Namespace, stack: ExitStack) -> TextIO:
    if args.infile:
        if args.gzip:
            return stack.enter_context(gzip.open(args.infile, "rt"))
        return stack.enter_context(open(args.infile))
    return sys.stdin


def open_output(args: argparse.Namespace, stack: ExitStack) -> TextIO:
    if args.outfile:
        return stack.enter_context(open(args.outfile, "w"))
    return sys.stdout
