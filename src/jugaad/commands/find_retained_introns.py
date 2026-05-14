"""Find transcript variants with retained introns in GFF3 or GTF annotations."""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict
from collections.abc import Sequence
from contextlib import ExitStack
from typing import TextIO

from jugaad.commands.common import add_gxf_arguments, currtime, open_input, open_output
from jugaad.parse_gxf_annots import (
    ExonsDict,
    Interval,
    get_gff3_exons,
    get_gtf_exons,
)

GeneKey = tuple[str, str, str]
TranscriptExons = dict[str, list[Interval]]
RetainedTranscript = tuple[str, str]

DESC_TEXT = (
    "This script parses input gff3/gtf file and reports transcript variants "
    "containing retained introns"
)


def add_arguments(parser: argparse.ArgumentParser) -> None:
    add_gxf_arguments(parser)


def build_parser(prog: str = "jugaad find_retained_introns") -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=prog, description=DESC_TEXT)
    add_arguments(parser)
    return parser


def register_parser(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "find_retained_introns",
        description=DESC_TEXT,
        help="find transcript variants with retained introns",
    )
    add_arguments(parser)
    parser.set_defaults(func=run)
    return parser


def group_exons_by_gene(exons_dict: ExonsDict) -> dict[GeneKey, TranscriptExons]:
    gene_exons: dict[GeneKey, TranscriptExons] = defaultdict(dict)
    for (chrom, strand, gene_id, transcript_id), exons in exons_dict.items():
        gene_exons[(chrom, strand, gene_id)][transcript_id] = sorted(exons)
    return gene_exons


def retained_exon_spans_reference_run(
    retained_exon: Interval,
    reference_exons: Sequence[Interval],
) -> bool:
    retained_start, retained_stop = retained_exon
    for start_index, first_exon in enumerate(reference_exons[:-1]):
        if first_exon[0] != retained_start:
            continue
        if first_exon[1] > retained_stop:
            continue

        for last_exon in reference_exons[start_index + 1 :]:
            if last_exon[0] > retained_stop:
                break
            if last_exon[1] > retained_stop:
                break
            if last_exon[1] == retained_stop:
                return True

    return False


def transcript_contains_retained_intron(
    candidate_exons: Sequence[Interval],
    reference_exons: Sequence[Interval],
) -> bool:
    if len(reference_exons) < 2:
        return False

    return any(
        retained_exon_spans_reference_run(candidate_exon, reference_exons)
        for candidate_exon in candidate_exons
    )


def find_retained_introns(exons_dict: ExonsDict) -> set[RetainedTranscript]:
    retained_transcripts: set[RetainedTranscript] = set()
    for (_chrom, _strand, gene_id), transcripts in group_exons_by_gene(
        exons_dict
    ).items():
        for candidate_id, candidate_exons in transcripts.items():
            for reference_id, reference_exons in transcripts.items():
                if candidate_id == reference_id:
                    continue
                if transcript_contains_retained_intron(
                    candidate_exons, reference_exons
                ):
                    retained_transcripts.add((candidate_id, gene_id))
                    break

    return retained_transcripts


def tabulate_retained_introns(
    retained_transcripts: set[RetainedTranscript],
    outfile: TextIO,
) -> None:
    tbl = csv.writer(outfile, delimiter="\t", lineterminator=os.linesep)
    tbl.writerow(["#transcript_id", "gene_id"])
    for transcript_id, gene_id in sorted(
        retained_transcripts,
        key=lambda transcript: (transcript[1], transcript[0]),
    ):
        tbl.writerow([transcript_id, gene_id])


def run(args: argparse.Namespace) -> int:
    tx_attrib = args.tx_attrib if args.tx_attrib else "transcript_id"
    gene_attrib = args.gene_attrib if args.gene_attrib else "GeneID"

    with ExitStack() as stack:
        infile = open_input(args, stack)
        outfile = open_output(args, stack)

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

        print(f"{currtime()} Finding retained introns...", file=sys.stderr)
        retained_transcripts = find_retained_introns(exons_dict)

        print(f"{currtime()} Writing output to file...", file=sys.stderr)
        tabulate_retained_introns(retained_transcripts, outfile)

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
