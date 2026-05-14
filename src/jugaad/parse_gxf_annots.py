"""Parse GFF3 and GTF annotation files."""

from __future__ import annotations

import csv
from collections import defaultdict
from typing import TextIO

csv.field_size_limit(10_000_000)

Interval = tuple[int, int]
TranscriptKey = tuple[str, str, str, str]
ExonsDict = defaultdict[TranscriptKey, list[Interval]]


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
