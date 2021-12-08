#!/usr/bin/env python3

import csv
import os
import sys
import re
import gzip
from datetime import datetime 
from collections import defaultdict
import argparse

# See http://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

def parse_args():
    desc_text = ("This script parses input gff3/gtf file and generates a table "
        "of introns")

    parser = argparse.ArgumentParser(description=desc_text)
    parser.add_argument('-i', '--infile', help="input file; default STDIN")
    parser.add_argument('-o', '--outfile', help="output file; default STDOUT")
    parser.add_argument('-z', '--gzip', help="input is gzip compressed", \
                        action='store_true')
    parser.add_argument('-f', '--format', help="input file type",
                        choices=['gtf', 'gff3'], required=True)
    parser.add_argument('--tx_attrib', help="attribute to use for transcript \
                        id; default is `transcript_id`", default='transcript_id')
    parser.add_argument('--gene_attrib', help="attribute to use for gene id; \
                        default is `GeneID`", default='GeneID')
    parser.add_argument('--splice_structures', help="write splice structures \
                        instead of introns to output file", action='store_true')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()

def currtime():
    return datetime.now().isoformat(timespec='seconds', sep=' ')


def get_gff3_exons(gff3_file, tx_attrib, gene_attrib):
    tbl = csv.reader(gff3_file, delimiter = '\t')
    exons_dict = defaultdict(list)
    for line in tbl:
        if not line[0].startswith('#'):
            assert len(line) == 9
            [ chrom, feat_source, feat_type,
             start, stop, score, strand,
             phase, attribs ] = line
            if feat_type == 'exon':
                attribs = process_attribs(attribs, '=')
                gene_id = attribs[gene_attrib]
                tx_id = attribs.get(tx_attrib, None)
                if not tx_id:
                    tx_id = attribs['Parent']
                start, stop = int(start), int(stop)
                exons_dict[(chrom, strand, gene_id, tx_id)].append((start, stop))
    gff3_file.close()
    return exons_dict

def get_gtf_exons(gtf_file, tx_attrib, gene_attrib):
    tbl = csv.reader(gtf_file, delimiter = '\t')
    exons_dict = defaultdict(list)
    for line in tbl:
        if not line[0].startswith('#'):
            assert len(line) == 9
            [chrom, feat_source, feat_type,
             start, stop, score, strand,
             phase, attribs ] = line
            if feat_type == 'exon':
                attribs = process_attribs(attribs, ' ')
                tx_id = attribs[tx_attrib]
                gene_id = attribs[gene_attrib]
                start, stop = map(int, [start, stop])
                exons_dict[(chrom, strand, gene_id, tx_id)].append((start, stop))
    gtf_file.close()
    return exons_dict

def process_attribs(attribs, delim):
    attrib_dict = {}
    attribs = list(filter(None, attribs.rstrip().split(';')))
    for attrib in attribs:
        k, v = attrib.lstrip().split(delim, 1)
        v = v.strip('"')
        if k == 'db_xref': ## gtf files have db_xref
            k, v = v.split(':', 1)
            attrib_dict[k] = v
        elif k == 'Dbxref': ## gff3 files have Dbxref
            xrefs = v.split(',')
            for xref in xrefs:
                k, v = xref.split(':', 1)
                attrib_dict[k] = v 
        elif k in attrib_dict:
            if isinstance(attrib_dict[k], str):
                attrib_dict[k] = [attrib_dict[k]]
            attrib_dict[k].append(v)
        else:
            attrib_dict[k] = v
    return attrib_dict

def get_introns(exons_dict):
    introns_dict = defaultdict(list)
    for tx_info, exons in exons_dict.items():
        if len(exons) > 1:
            [ chrom, strand, gene_id, transcript_id ] = tx_info
            exons.sort()
            i = 0
            while (i + 1) < len(exons):
                introns = (exons[i][1] + 1, exons[i+1][0] - 1)
                i = i + 1
                introns_dict[(chrom, strand, gene_id, transcript_id)].append(introns)
    return introns_dict

def tabulate_introns(introns_dict, outfile):
    ## introns in the table are ordered based on their position
    ## in the transcript so the order of introns is dependent
    ## on the strand the transcript is on
    tbl = csv.writer(outfile,
                     delimiter = '\t',
                     lineterminator = os.linesep)
    tbl.writerow(['#intron_id', 'gene_id', 'refseq_acc', 
        'intron_num', 'intron_ct'])
    for rs_info, introns in introns_dict.items():
        [
            chrom,
            strand,
            gene_id,
            rs_acc
        ] = rs_info
        if strand == '+':
            introns = sorted(introns)
        elif strand == '-':
            introns = sorted(introns, reverse = True)
        num_introns = len(introns)
        for intron in introns:
            intron_id = '|'.join([chrom, str(intron[0]), str(intron[1]), strand])
            tbl.writerow([intron_id,
                         gene_id,
                         rs_acc,
                         introns.index(intron) + 1,
                         num_introns])
    outfile.close()

def tabulate_splice_structures(introns_dict, outfile):
    ## introns in a splice structure always go left to right
    ## irrespective of the strand transcript is on
    tbl = csv.writer(outfile,
                     delimiter = '\t',
                     lineterminator = os.linesep)
    tbl.writerow(['#transcript_id', 'gene_id', 'splice_structure'])
    for tx_info, introns in introns_dict.items():
        [ chrom, strand, gene_id, transcript_id ] = tx_info
        splice_structure = []
        for intron in introns:
            intron_id = '|'.join([chrom, str(intron[0]), str(intron[1]), strand])
            splice_structure.append(intron_id)
        tbl.writerow([transcript_id, gene_id, ','.join(splice_structure)])
    outfile.close()


## setup args
args = parse_args()

if args.infile:
    if args.gzip:
        infile = gzip.open(args.infile, 'rt')
    else:
        infile = open(args.infile, 'r')
else:
    infile = sys.stdin

if args.outfile:
    outfile = open(args.outfile, 'w')
else:
    outfile = sys.stdout

tx_attrib = args.tx_attrib if args.tx_attrib else 'transcript_id'
gene_attrib = args.gene_attrib if args.gene_attrib else 'GeneID'

## create exon dictionary
print('{} Parsing input file and creating exon dictionary...' 
    .format(currtime()), file=sys.stderr)
if args.format == 'gff3':
    exons_dict = get_gff3_exons(infile, tx_attrib, gene_attrib)
elif args.format == 'gtf':
    exons_dict = get_gtf_exons(infile, tx_attrib, gene_attrib)

## create intron dictionary
print('{} Computing introns...' 
    .format(currtime()), file=sys.stderr)
introns_dict = get_introns(exons_dict)

## write output to file
print('{} Writing output to file...' 
    .format(currtime()), file=sys.stderr)
if args.splice_structures:
    tabulate_splice_structures(introns_dict, outfile)
else:
    tabulate_introns(introns_dict, outfile)
