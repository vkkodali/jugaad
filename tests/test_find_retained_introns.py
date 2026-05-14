from __future__ import annotations

import gzip
import io
import sys

import pytest

from jugaad.cli import main as cli_main


GFF3_RETAINED_DATA = """##gff-version 3
chr1\tRefSeq\texon\t100\t200\t.\t+\t.\tID=exon1;Parent=TX1;Dbxref=GeneID:GENE1
chr1\tRefSeq\texon\t300\t400\t.\t+\t.\tID=exon2;Parent=TX1;Dbxref=GeneID:GENE1
chr1\tRefSeq\texon\t500\t600\t.\t+\t.\tID=exon3;Parent=TX1;Dbxref=GeneID:GENE1
chr1\tRefSeq\texon\t800\t900\t.\t+\t.\tID=exon4;Parent=TX1;Dbxref=GeneID:GENE1
chr1\tRefSeq\texon\t100\t200\t.\t+\t.\tID=exon5;Parent=TX2;Dbxref=GeneID:GENE1
chr1\tRefSeq\texon\t300\t600\t.\t+\t.\tID=exon6;Parent=TX2;Dbxref=GeneID:GENE1
chr1\tRefSeq\texon\t800\t900\t.\t+\t.\tID=exon7;Parent=TX2;Dbxref=GeneID:GENE1
"""

GFF3_MULTI_SPAN_DATA = """##gff-version 3
chr1\tRefSeq\texon\t100\t150\t.\t+\t.\tID=exon1;transcript_id=TX1;GeneID=GENE1
chr1\tRefSeq\texon\t200\t250\t.\t+\t.\tID=exon2;transcript_id=TX1;GeneID=GENE1
chr1\tRefSeq\texon\t300\t350\t.\t+\t.\tID=exon3;transcript_id=TX1;GeneID=GENE1
chr1\tRefSeq\texon\t400\t450\t.\t+\t.\tID=exon4;transcript_id=TX1;GeneID=GENE1
chr1\tRefSeq\texon\t500\t550\t.\t+\t.\tID=exon5;transcript_id=TX1;GeneID=GENE1
chr1\tRefSeq\texon\t100\t150\t.\t+\t.\tID=exon6;transcript_id=TX3;GeneID=GENE1
chr1\tRefSeq\texon\t200\t250\t.\t+\t.\tID=exon7;transcript_id=TX3;GeneID=GENE1
chr1\tRefSeq\texon\t300\t350\t.\t+\t.\tID=exon8;transcript_id=TX3;GeneID=GENE1
chr1\tRefSeq\texon\t400\t450\t.\t+\t.\tID=exon9;transcript_id=TX3;GeneID=GENE1
chr1\tRefSeq\texon\t500\t550\t.\t+\t.\tID=exon10;transcript_id=TX3;GeneID=GENE1
chr1\tRefSeq\texon\t100\t150\t.\t+\t.\tID=exon11;transcript_id=TX2;GeneID=GENE1
chr1\tRefSeq\texon\t200\t450\t.\t+\t.\tID=exon12;transcript_id=TX2;GeneID=GENE1
chr1\tRefSeq\texon\t500\t550\t.\t+\t.\tID=exon13;transcript_id=TX2;GeneID=GENE1
"""

GFF3_NON_MATCH_DATA = """##gff-version 3
chr1\tRefSeq\texon\t100\t200\t.\t+\t.\tID=exon1;transcript_id=ONE_REF;GeneID=GENE1
chr1\tRefSeq\texon\t300\t400\t.\t+\t.\tID=exon2;transcript_id=ONE_REF;GeneID=GENE1
chr1\tRefSeq\texon\t100\t200\t.\t+\t.\tID=exon3;transcript_id=ONE_CAND;GeneID=GENE1
chr2\tRefSeq\texon\t1000\t1100\t.\t+\t.\tID=exon4;transcript_id=BAD_REF;GeneID=GENE2
chr2\tRefSeq\texon\t1200\t1300\t.\t+\t.\tID=exon5;transcript_id=BAD_REF;GeneID=GENE2
chr2\tRefSeq\texon\t1000\t1290\t.\t+\t.\tID=exon6;transcript_id=BAD_CAND;GeneID=GENE2
chr3\tRefSeq\texon\t100\t150\t.\t+\t.\tID=exon7;transcript_id=GENE_REF;GeneID=GENE3A
chr3\tRefSeq\texon\t200\t250\t.\t+\t.\tID=exon8;transcript_id=GENE_REF;GeneID=GENE3A
chr3\tRefSeq\texon\t100\t250\t.\t+\t.\tID=exon9;transcript_id=GENE_CAND;GeneID=GENE3B
chr4\tRefSeq\texon\t100\t150\t.\t+\t.\tID=exon10;transcript_id=CHROM_REF;GeneID=GENE4
chr4\tRefSeq\texon\t200\t250\t.\t+\t.\tID=exon11;transcript_id=CHROM_REF;GeneID=GENE4
chr5\tRefSeq\texon\t100\t250\t.\t+\t.\tID=exon12;transcript_id=CHROM_CAND;GeneID=GENE4
chr6\tRefSeq\texon\t100\t150\t.\t+\t.\tID=exon13;transcript_id=STRAND_REF;GeneID=GENE5
chr6\tRefSeq\texon\t200\t250\t.\t+\t.\tID=exon14;transcript_id=STRAND_REF;GeneID=GENE5
chr6\tRefSeq\texon\t100\t250\t.\t-\t.\tID=exon15;transcript_id=STRAND_CAND;GeneID=GENE5
"""


def run_cli(monkeypatch: pytest.MonkeyPatch, argv: list[str], stdin: str = ""):
    stdout = io.StringIO()
    stderr = io.StringIO()
    monkeypatch.setattr(sys, "stdin", io.StringIO(stdin))
    monkeypatch.setattr(sys, "stdout", stdout)
    monkeypatch.setattr(sys, "stderr", stderr)

    exit_code = cli_main(argv)
    return exit_code, stdout.getvalue(), stderr.getvalue()


def test_gff3_reports_retained_transcript_from_stdin(monkeypatch):
    exit_code, stdout, stderr = run_cli(
        monkeypatch,
        ["find_retained_introns", "-f", "gff3"],
        GFF3_RETAINED_DATA,
    )

    assert exit_code == 0
    assert stdout == "#transcript_id\tgene_id\nTX2\tGENE1\n"
    assert "Parsing input file and creating exon dictionary" in stderr
    assert "Finding retained introns" in stderr
    assert "Writing output to file" in stderr


def test_multi_exon_span_reports_transcript_once(monkeypatch):
    exit_code, stdout, _stderr = run_cli(
        monkeypatch,
        ["find_retained_introns", "-f", "gff3"],
        GFF3_MULTI_SPAN_DATA,
    )

    assert exit_code == 0
    assert stdout == "#transcript_id\tgene_id\nTX2\tGENE1\n"


def test_non_matches_are_ignored(monkeypatch):
    exit_code, stdout, _stderr = run_cli(
        monkeypatch,
        ["find_retained_introns", "-f", "gff3"],
        GFF3_NON_MATCH_DATA,
    )

    assert exit_code == 0
    assert stdout == "#transcript_id\tgene_id\n"


def test_gzip_input(monkeypatch, tmp_path):
    infile = tmp_path / "input.gff3.gz"
    with gzip.open(infile, "wt") as gff3_file:
        gff3_file.write(GFF3_RETAINED_DATA)

    exit_code, stdout, _stderr = run_cli(
        monkeypatch,
        ["find_retained_introns", "-f", "gff3", "-z", "-i", str(infile)],
    )

    assert exit_code == 0
    assert stdout == "#transcript_id\tgene_id\nTX2\tGENE1\n"


def test_outfile_and_attribute_overrides(monkeypatch, tmp_path):
    infile = tmp_path / "input.gff3"
    outfile = tmp_path / "retained_introns.tsv"
    infile.write_text(
        "chr4\tRefSeq\texon\t10\t20\t.\t+\t.\tID=exon1;custom_tx=TX_REF;"
        "custom_gene=GENE4\n"
        "chr4\tRefSeq\texon\t30\t40\t.\t+\t.\tID=exon2;custom_tx=TX_REF;"
        "custom_gene=GENE4\n"
        "chr4\tRefSeq\texon\t10\t40\t.\t+\t.\tID=exon3;custom_tx=TX_RET;"
        "custom_gene=GENE4\n"
    )

    exit_code, stdout, _stderr = run_cli(
        monkeypatch,
        [
            "find_retained_introns",
            "-f",
            "gff3",
            "-i",
            str(infile),
            "-o",
            str(outfile),
            "--tx_attrib",
            "custom_tx",
            "--gene_attrib",
            "custom_gene",
        ],
    )

    assert exit_code == 0
    assert stdout == ""
    assert outfile.read_text() == "#transcript_id\tgene_id\nTX_RET\tGENE4\n"


def test_gtf_input_uses_geneid_default(monkeypatch):
    gtf_data = (
        'chr2\tRefSeq\texon\t100\t150\t.\t+\t.\tgene_id "ignored"; '
        'transcript_id "TX1"; db_xref "GeneID:GENE2";\n'
        'chr2\tRefSeq\texon\t250\t300\t.\t+\t.\tgene_id "ignored"; '
        'transcript_id "TX1"; db_xref "GeneID:GENE2";\n'
        'chr2\tRefSeq\texon\t100\t300\t.\t+\t.\tgene_id "ignored"; '
        'transcript_id "TX2"; db_xref "GeneID:GENE2";\n'
    )

    exit_code, stdout, _stderr = run_cli(
        monkeypatch,
        ["find_retained_introns", "-f", "gtf"],
        gtf_data,
    )

    assert exit_code == 0
    assert stdout == "#transcript_id\tgene_id\nTX2\tGENE2\n"


def test_subcommand_help(monkeypatch):
    stdout = io.StringIO()
    stderr = io.StringIO()
    monkeypatch.setattr(sys, "stdout", stdout)
    monkeypatch.setattr(sys, "stderr", stderr)

    with pytest.raises(SystemExit) as excinfo:
        cli_main(["find_retained_introns", "--help"])

    assert excinfo.value.code == 0
    assert "usage: jugaad find_retained_introns" in stdout.getvalue()
    assert "input file type" in stdout.getvalue()
    assert "retained introns" in stdout.getvalue()
    assert stderr.getvalue() == ""
