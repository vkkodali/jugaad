"""Generate intron tables from GFF3 or GTF annotations."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
from collections import defaultdict
from collections.abc import Sequence
from contextlib import ExitStack
from datetime import datetime
from typing import TextIO

csv.field_size_limit(10_000_000)

Interval = tuple[int, int]
TranscriptKey = tuple[str, str, str, str]
ExonsDict = defaultdict[TranscriptKey, list[Interval]]
IntronsDict = defaultdict[TranscriptKey, list[Interval]]

DESC_TEXT = "This script parses input gff3/gtf file and generates a table of introns"


def currtime() -> str:
    return datetime.now().isoformat(timespec="seconds", sep=" ")


def add_arguments(parser: argparse.ArgumentParser) -> None:
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
    parser.add_argument(
        "--splice_structures",
        help="write splice structures instead of introns to output file",
        action="store_true",
    )


def build_parser(prog: str = "jugaad gff3_to_introns") -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=prog, description=DESC_TEXT)
    add_arguments(parser)
    return parser


def register_parser(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "gff3_to_introns",
        description=DESC_TEXT,
        help="generate a table of introns from a GFF3/GTF file",
    )
    add_arguments(parser)
    parser.set_defaults(func=run)
    return parser


def process_attribs(attribs: str, delim: str) -> dict[str, str | list[str]]:
    attrib_dict: dict[str, str | list[str]] = {}
    attrib_parts = list(filter(None, attribs.rstrip().split(";")))
    for attrib in attrib_parts:
        key, value = attrib.lstrip().split(delim, 1)
        value = value.strip('"')
        if key == "db_xref":
            key, value = value.split(":", 1)
            attrib_dict[key] = value
        elif key == "Dbxref":
            xrefs = value.split(",")
            for xref in xrefs:
                key, value = xref.split(":", 1)
                attrib_dict[key] = value
        elif key in attrib_dict:
            if isinstance(attrib_dict[key], str):
                attrib_dict[key] = [attrib_dict[key]]
            attrib_dict[key].append(value)
        else:
            attrib_dict[key] = value
    return attrib_dict


def get_gff3_exons(
    gff3_file: TextIO,
    tx_attrib: str,
    gene_attrib: str,
) -> ExonsDict:
    tbl = csv.reader(gff3_file, delimiter="\t")
    exons_dict: ExonsDict = defaultdict(list)
    for line in tbl:
        if not line or line[0].startswith("#"):
            continue
        assert len(line) == 9
        [
            chrom,
            _feat_source,
            feat_type,
            start,
            stop,
            _score,
            strand,
            _phase,
            attribs,
        ] = line
        if feat_type == "exon":
            attrib_dict = process_attribs(attribs, "=")
            gene_id = attrib_dict[gene_attrib]
            tx_id = attrib_dict.get(tx_attrib)
            if not tx_id:
                tx_id = attrib_dict["Parent"]
            start_int, stop_int = int(start), int(stop)
            exons_dict[(chrom, strand, gene_id, tx_id)].append((start_int, stop_int))
    return exons_dict


def get_gtf_exons(
    gtf_file: TextIO,
    tx_attrib: str,
    gene_attrib: str,
) -> ExonsDict:
    tbl = csv.reader(gtf_file, delimiter="\t")
    exons_dict: ExonsDict = defaultdict(list)
    for line in tbl:
        if not line or line[0].startswith("#"):
            continue
        assert len(line) == 9
        [
            chrom,
            _feat_source,
            feat_type,
            start,
            stop,
            _score,
            strand,
            _phase,
            attribs,
        ] = line
        if feat_type == "exon":
            attrib_dict = process_attribs(attribs, " ")
            tx_id = attrib_dict[tx_attrib]
            gene_id = attrib_dict[gene_attrib]
            start_int, stop_int = map(int, [start, stop])
            exons_dict[(chrom, strand, gene_id, tx_id)].append((start_int, stop_int))
    return exons_dict


def get_introns(exons_dict: ExonsDict) -> IntronsDict:
    introns_dict: IntronsDict = defaultdict(list)
    for tx_info, exons in exons_dict.items():
        if len(exons) > 1:
            chrom, strand, gene_id, transcript_id = tx_info
            exons.sort()
            index = 0
            while (index + 1) < len(exons):
                intron = (exons[index][1] + 1, exons[index + 1][0] - 1)
                index += 1
                introns_dict[(chrom, strand, gene_id, transcript_id)].append(intron)
    return introns_dict


def tabulate_introns(introns_dict: IntronsDict, outfile: TextIO) -> None:
    tbl = csv.writer(outfile, delimiter="\t", lineterminator=os.linesep)
    tbl.writerow(["#intron_id", "gene_id", "refseq_acc", "intron_num", "intron_ct"])
    for rs_info, introns in introns_dict.items():
        chrom, strand, gene_id, rs_acc = rs_info
        if strand == "+":
            ordered_introns = sorted(introns)
        elif strand == "-":
            ordered_introns = sorted(introns, reverse=True)
        else:
            ordered_introns = introns

        num_introns = len(ordered_introns)
        for intron_num, intron in enumerate(ordered_introns, start=1):
            intron_id = "|".join([chrom, str(intron[0]), str(intron[1]), strand])
            tbl.writerow([intron_id, gene_id, rs_acc, intron_num, num_introns])


def tabulate_splice_structures(introns_dict: IntronsDict, outfile: TextIO) -> None:
    tbl = csv.writer(outfile, delimiter="\t", lineterminator=os.linesep)
    tbl.writerow(["#transcript_id", "gene_id", "splice_structure"])
    for tx_info, introns in introns_dict.items():
        chrom, strand, gene_id, transcript_id = tx_info
        splice_structure = []
        for intron in introns:
            intron_id = "|".join([chrom, str(intron[0]), str(intron[1]), strand])
            splice_structure.append(intron_id)
        tbl.writerow([transcript_id, gene_id, ",".join(splice_structure)])


def _open_input(args: argparse.Namespace, stack: ExitStack) -> TextIO:
    if args.infile:
        if args.gzip:
            return stack.enter_context(gzip.open(args.infile, "rt"))
        return stack.enter_context(open(args.infile))
    return sys.stdin


def _open_output(args: argparse.Namespace, stack: ExitStack) -> TextIO:
    if args.outfile:
        return stack.enter_context(open(args.outfile, "w"))
    return sys.stdout


def run(args: argparse.Namespace) -> int:
    tx_attrib = args.tx_attrib if args.tx_attrib else "transcript_id"
    gene_attrib = args.gene_attrib if args.gene_attrib else "GeneID"

    with ExitStack() as stack:
        infile = _open_input(args, stack)
        outfile = _open_output(args, stack)

        print(
            f"{currtime()} Parsing input file and creating exon dictionary...",
            file=sys.stderr,
        )
        if args.format == "gff3":
            exons_dict = get_gff3_exons(infile, tx_attrib, gene_attrib)
        elif args.format == "gtf":
            exons_dict = get_gtf_exons(infile, tx_attrib, gene_attrib)
        else:
            raise ValueError(f"Unsupported format: {args.format}")

        print(f"{currtime()} Computing introns...", file=sys.stderr)
        introns_dict = get_introns(exons_dict)

        print(f"{currtime()} Writing output to file...", file=sys.stderr)
        if args.splice_structures:
            tabulate_splice_structures(introns_dict, outfile)
        else:
            tabulate_introns(introns_dict, outfile)

    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args_list = list(sys.argv[1:] if argv is None else argv)
    if not args_list:
        parser.print_help(sys.stderr)
        return 1
    return run(parser.parse_args(args_list))


if __name__ == "__main__":
    raise SystemExit(main())
