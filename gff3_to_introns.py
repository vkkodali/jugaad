#!/usr/bin/env python

import csv
import os
import sys
import gzip
from collections import defaultdict
import argparse

# See http://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

parser = argparse.ArgumentParser(description ="""This script parses input
                gff3 file and generates a table of introns""")
parser.add_argument('-i', '--infile', help="input file; default STDIN")
parser.add_argument('-o', '--outfile', help="output file; default STDOUT")
parser.add_argument('-z', '--gzip', help="input is gzip compressed", \
                    action = 'store_true')
args = parser.parse_args()

if args.infile:
    if args.gzip:
        gff3_file = gzip.open(args.infile, 'rt')
    else:
        gff3_file = open(args.infile, 'r')
else:
    gff3_file = sys.stdin

if args.outfile:
    introns_file = open(args.outfile, 'w')
else:
    introns_file = sys.stdout

def get_exons(gff3_file):
    tbl = csv.reader(gff3_file, delimiter = '\t')
    exons_dict = defaultdict(list)
    for line in tbl:
        if not line[0].startswith('#'):
            [
                chrom,
                feat_source,
                feat_type,
                start,
                stop,
                score,
                strand,
                phase,
                attribs
            ] = line
            if (feat_type == "exon"
                and 'transcript_id=' in attribs):
                new_attribs = process_attribs(attribs)
                if 'GeneID' in new_attribs:
                    gene_id = new_attribs['GeneID']
                else:
                    gene_id = new_attribs['gene']
                tx = new_attribs['transcript_id']
                start, stop = int(start), int(stop)
                exons_dict[(chrom, strand, gene_id, tx)].append((start, stop))
            elif(feat_type == "exon"
                 and 'Parent=transcript' in attribs):
                attribs = attribs.split(';')
                new_attribs = {}
                for attrib in attribs:
                    attrib = attrib.split('=')
                    new_attribs[attrib[0]] = attrib[1]
                gene_id = new_attribs['Parent'].strip('transcript:')
                tx = new_attribs['Parent'].strip('transcript:')
                start, stop = int(start), int(stop)
                exons_dict[(chrom, strand, gene_id, tx)].append((start, stop))
    gff3_file.close()
    return exons_dict

def process_attribs(attribs):
    new_attribs = {}
    attribs = list(filter(None, attribs.split(';'))) ## removes empty strings, needed because some gff3 lines have ";;"
    for attrib in attribs:
        k, v = attrib.split('=')
        if k == 'Dbxref':
            xrefs = v.split(',')
            for xref in xrefs:
                terms = xref.split(':')
                new_attribs[terms[-2]] = terms[-1]
        else:
            new_attribs[k] = v
    return new_attribs

def get_introns(exons_dict):
    introns_dict = defaultdict(set)
    for rs_info, exons in exons_dict.items():
        if len(exons) > 1:
            [
                chrom,
                strand,
                gene_id,
                rs
            ] = rs_info
            exons = sorted(exons)
            i = 0
            while (i + 1) < len(exons):
                introns = (exons[i][1] + 1, exons[i+1][0] - 1)
                i = i + 1
                introns_dict[(chrom, strand, gene_id, rs)].add(introns)
    return introns_dict

def tabulate_introns(introns_dict, introns_file):
    tbl = csv.writer(introns_file,
                     delimiter = '\t',
                     lineterminator = os.linesep)
    tbl.writerow(['#intron_id',
                  'gene_id',
                  'refseq_acc',
                  'intron_num',
                  'intron_ct'])
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
    introns_file.close()

exons_dict = get_exons(gff3_file)
introns_dict = get_introns(exons_dict)
tabulate_introns(introns_dict, introns_file)
