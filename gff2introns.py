#!/usr/bin/env python

import csv
import argparse
from collections import defaultdict

def getintrons(infile):
    with open(infile, 'r') as f:
        tbl = csv.reader(f, delimiter = '\t')
        intron_dict = defaultdict(list)
        old_id = 'None'
        for line in tbl:
            if not line[0].startswith('#'):
                [   chrom,
                    pa,
                    feat_type,
                    start,
                    stop,
                    score,
                    strand,
                    phase,
                    attribs
                ] = line
                [   gff_id,
                    dbxref,
                    pa_acc,
                    seq,
                    n_obs,
                    n_sams,
                    n_maps,
                    n_locs
                ] = attribs.split(';')
                if feat_type == 'CDS':
                    if gff_id != old_id:
                        old_id, old_start, old_stop = gff_id, start, stop
                    elif gff_id == old_id:
                        if old_stop < start:
                            introns = ((int(old_stop) + 1), (int(start) - 1))
                        elif old_stop > start:
                            introns = ((int(stop) + 1), (int(old_start) - 1))
                        intron_dict[introns].append([chrom, strand, n_obs])
    return intron_dict

def writegff(intron_dict):
    with open(outfile, 'w') as f:
        tbl = csv.writer(f, delimiter = '\t')
        gff_id = 0
        for introns in intron_dict:
            start, stop = introns
            count = 0
            for item in intron_dict[introns]:
                chrom, strand, n_obs = item
                count = count + int(n_obs.split("=")[1])
            attribs = [ "ID=id"+str(gff_id),
                        "count="+str(count)
                      ]
            gff3_line = [   chrom,
                            'PeptideAtlas',
                            'intron',
                            start,
                            stop,
                            count,
                            strand,
                            '0',
                            ";".join(map(str, attribs))
                        ]
            tbl.writerow(gff3_line)
            gff_id = gff_id + 1

parser = argparse.ArgumentParser(description='''This script parses a GFF3
file, extracts introns and generates an output file in GFF3 format''')
parser.add_argument('infile', help='input GFF3 file')
parser.add_argument('outfile', help='output GFF3 file')
args = parser.parse_args()

infile = args.infile
outfile = args.outfile

intron_dict = getintrons(infile)
writegff(intron_dict)
