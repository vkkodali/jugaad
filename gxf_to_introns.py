#!/usr/bin/env python3

import argparse
import os
import csv
import sys
from datetime import datetime
from collections import defaultdict

try:
    import gffutils
except ModuleNotFoundError as e :
    print(e, 'Please install gffutils', sep = '\n')
    sys.exit(1)

def currtime():
    return datetime.now().isoformat(timespec='seconds', sep=' ')

def parse_args():
    desc_text = ("This script parses input gff3/gtf file to produce splice "
        "structure table")
    parser = argparse.ArgumentParser(description = desc_text)
    parser.add_argument('-i', '--infile', help="input file; can be gzipped",
        required = True, )
    parser.add_argument('-o', '--outfile', help="output file; default STDOUT")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()

def create_db(gxfin):
    dbfilename = '/tmp/' + datetime.now().strftime('%Y%m%d%H%M%S%f') + '_gffutils.db'
    db = gffutils.create_db(gxfin, 
                dbfn=dbfilename, 
                force=True,
                disable_infer_genes=True,
                disable_infer_transcripts=True) 
    return db

def add_introns(db):
    db.update(db.create_introns())

def create_intron_dict(db):
    d = defaultdict(list)
    for intron in db.features_of_type('intron'):
        intron_id = [intron.chrom, intron.start, intron.stop, intron.strand]
        txid = intron.attributes['transcript_id'][0]
        gxid = intron.attributes['gene_id'][0]
        d[(gxid, txid)].append(intron_id)
    return d

def tabulate_splice_structures(fo, d):
    tbl = csv.writer(fo, delimiter = '\t', lineterminator = os.linesep)
    tbl.writerow(['#gene_id', 'transcript_id', 'intron_id'])

    for txinfo, introns in d.items():
        txinfo = list(txinfo)
        introns.sort(key = lambda x: x[1]) # https://stackoverflow.com/a/4174956
        for i in range(0, len(introns)):
            for j in range(0, len(introns[i])):
                if isinstance(introns[i][j], int):
                    introns[i][j] = str(introns[i][j])
        introns = ['|'.join(i) for i in introns]
        tbl.writerow(txinfo + [','.join(introns)])

## parse arguments
args = parse_args()
gxfin = args.infile
if args.outfile:
    fo = open(args.outfile, 'wt')
else:
    fo = sys.stdout

# create db
print('{} Creating database from the input file...' 
    .format(currtime()), file=sys.stderr)
db = create_db(gxfin)

# compute intron features and update db
print('{} Computing introns...' 
    .format(currtime()), file=sys.stderr)
add_introns(db)

# create a dict with introns for each tx
print('{} Building splice structures...' 
    .format(currtime()), file=sys.stderr)
int_feats = create_intron_dict(db)

# compute splice structures and write to file
print('{} Writing output to file...' 
    .format(currtime()), file=sys.stderr)
tabulate_splice_structures(fo, int_feats)

