from __future__ import annotations

import gzip
import io
import sys

import pytest

from jugaad.cli import main as cli_main
from jugaad.commands.gff3_to_introns import process_attribs


GFF3_PARENT_DATA = """##gff-version 3
chr1\tRefSeq\texon\t100\t200\t.\t+\t.\tID=exon1;Parent=TXP1;Dbxref=GeneID:GENE1
chr1\tRefSeq\texon\t300\t400\t.\t+\t.\tID=exon2;Parent=TXP1;Dbxref=GeneID:GENE1
"""

GFF3_MINUS_DATA = """##gff-version 3
chr3\tRefSeq\texon\t100\t150\t.\t-\t.\tID=exon1;transcript_id=TX3;GeneID=GENE3
chr3\tRefSeq\texon\t200\t250\t.\t-\t.\tID=exon2;transcript_id=TX3;GeneID=GENE3
chr3\tRefSeq\texon\t300\t350\t.\t-\t.\tID=exon3;transcript_id=TX3;GeneID=GENE3
"""


def run_cli(monkeypatch: pytest.MonkeyPatch, argv: list[str], stdin: str = ""):
    stdout = io.StringIO()
    stderr = io.StringIO()
    monkeypatch.setattr(sys, "stdin", io.StringIO(stdin))
    monkeypatch.setattr(sys, "stdout", stdout)
    monkeypatch.setattr(sys, "stderr", stderr)

    exit_code = cli_main(argv)
    return exit_code, stdout.getvalue(), stderr.getvalue()


def test_gff3_stdin_stdout_parent_fallback_and_dbxref(monkeypatch):
    exit_code, stdout, stderr = run_cli(
        monkeypatch,
        ["gff3_to_introns", "-f", "gff3"],
        GFF3_PARENT_DATA,
    )

    assert exit_code == 0
    assert stdout == (
        "#intron_id\tgene_id\trefseq_acc\tintron_num\tintron_ct\n"
        "chr1|201|299|+\tGENE1\tTXP1\t1\t1\n"
    )
    assert "Parsing input file and creating exon dictionary" in stderr
    assert "Computing introns" in stderr
    assert "Writing output to file" in stderr


def test_gtf_input_uses_db_xref_and_duplicate_attributes(monkeypatch):
    attrs = process_attribs(
        'gene_id "GENE2"; transcript_id "TX2"; tag "one"; tag "two"; '
        'db_xref "GeneID:GENE2";',
        " ",
    )
    assert attrs["tag"] == ["one", "two"]
    assert attrs["GeneID"] == "GENE2"

    gtf_data = (
        'chr2\tRefSeq\texon\t100\t150\t.\t+\t.\tgene_id "ignored"; '
        'transcript_id "TX2"; db_xref "GeneID:GENE2";\n'
        'chr2\tRefSeq\texon\t250\t300\t.\t+\t.\tgene_id "ignored"; '
        'transcript_id "TX2"; db_xref "GeneID:GENE2";\n'
    )
    exit_code, stdout, _stderr = run_cli(
        monkeypatch,
        ["gff3_to_introns", "-f", "gtf"],
        gtf_data,
    )

    assert exit_code == 0
    assert stdout == (
        "#intron_id\tgene_id\trefseq_acc\tintron_num\tintron_ct\n"
        "chr2|151|249|+\tGENE2\tTX2\t1\t1\n"
    )


def test_gzip_input(monkeypatch, tmp_path):
    infile = tmp_path / "input.gff3.gz"
    with gzip.open(infile, "wt") as gff3_file:
        gff3_file.write(GFF3_PARENT_DATA)

    exit_code, stdout, _stderr = run_cli(
        monkeypatch,
        ["gff3_to_introns", "-f", "gff3", "-z", "-i", str(infile)],
    )

    assert exit_code == 0
    assert "chr1|201|299|+\tGENE1\tTXP1\t1\t1\n" in stdout


def test_outfile_and_attribute_overrides(monkeypatch, tmp_path):
    infile = tmp_path / "input.gff3"
    outfile = tmp_path / "introns.tsv"
    infile.write_text(
        "chr4\tRefSeq\texon\t10\t20\t.\t+\t.\tID=exon1;custom_tx=TX4;"
        "custom_gene=GENE4\n"
        "chr4\tRefSeq\texon\t30\t40\t.\t+\t.\tID=exon2;custom_tx=TX4;"
        "custom_gene=GENE4\n"
    )

    exit_code, stdout, _stderr = run_cli(
        monkeypatch,
        [
            "gff3_to_introns",
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
    assert outfile.read_text() == (
        "#intron_id\tgene_id\trefseq_acc\tintron_num\tintron_ct\n"
        "chr4|21|29|+\tGENE4\tTX4\t1\t1\n"
    )


def test_minus_strand_intron_numbering(monkeypatch):
    exit_code, stdout, _stderr = run_cli(
        monkeypatch,
        ["gff3_to_introns", "-f", "gff3"],
        GFF3_MINUS_DATA,
    )

    assert exit_code == 0
    assert stdout == (
        "#intron_id\tgene_id\trefseq_acc\tintron_num\tintron_ct\n"
        "chr3|251|299|-\tGENE3\tTX3\t1\t2\n"
        "chr3|151|199|-\tGENE3\tTX3\t2\t2\n"
    )


def test_splice_structures_output_left_to_right(monkeypatch):
    exit_code, stdout, _stderr = run_cli(
        monkeypatch,
        ["gff3_to_introns", "-f", "gff3", "--splice_structures"],
        GFF3_MINUS_DATA,
    )

    assert exit_code == 0
    assert stdout == (
        "#transcript_id\tgene_id\tsplice_structure\n"
        "TX3\tGENE3\tchr3|151|199|-,chr3|251|299|-\n"
    )


def test_no_subcommand_prints_help_to_stderr(monkeypatch):
    exit_code, stdout, stderr = run_cli(monkeypatch, [])

    assert exit_code == 1
    assert stdout == ""
    assert "usage: jugaad" in stderr
    assert "gff3_to_introns" in stderr


def test_subcommand_help(monkeypatch):
    stdout = io.StringIO()
    stderr = io.StringIO()
    monkeypatch.setattr(sys, "stdout", stdout)
    monkeypatch.setattr(sys, "stderr", stderr)

    with pytest.raises(SystemExit) as excinfo:
        cli_main(["gff3_to_introns", "--help"])

    assert excinfo.value.code == 0
    assert "usage: jugaad gff3_to_introns" in stdout.getvalue()
    assert "input file type" in stdout.getvalue()
    assert stderr.getvalue() == ""
