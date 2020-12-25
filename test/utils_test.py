from __future__ import print_function

import os.path as op

import clodius.utils as cu
from clodius.tiles.utils import abs2genome_fn


def test_infer_filetype():
    assert cu.infer_filetype("blah.gff") == "gff"
    assert cu.infer_filetype("blah.gff.gz") == "gff"
    assert cu.infer_filetype("blah.xyz") == None
    assert cu.infer_filetype("blah.bam") == "bam"
    assert cu.infer_filetype("blah.bed.bgz") == "bedfile"
    assert cu.infer_filetype("blah.bed") == "bedfile"


def test_infer_datatype():
    assert cu.infer_datatype("gff") == "bedlike"
    assert cu.infer_datatype("cooler") == "matrix"
    assert cu.infer_datatype("bedfile") == "bedlike"
    assert cu.infer_datatype("bam") == "reads"


def test_abs2genome_fn():
    fai_filename = op.join(
        "data", "GCA_000350705.1_Esch_coli_KTE11_V1_genomic.short.fna.fai"
    )
    sections = list(abs2genome_fn(fai_filename, 0, 1000))

    assert len(sections) == 3
    assert sections[0].end == 640
