import clodius.chromosomes as cs
import clodius.tiles.bedpe as ctbp
import os.path as op

testdir = op.realpath(op.dirname(__file__))

import numpy as np
import pytest


@pytest.mark.parametrize(
    "filename,header",
    [
        (
            "isidro.bedpe",
            "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsv_id\tpe_support\tstrand1\tstrand2\tsvclass\tsvmethod",
        ),
        ("hg19_myc.bedpe", ""),
    ],
)
def test_bedpe_tileset_info(filename, header):
    input_file = op.join(testdir, "sample_data", filename)
    chromsizes_fn = op.join(testdir, "sample_data", "b37.chrom.sizes")

    chromsizes = cs.chromsizes_as_series(chromsizes_fn)
    tileset_info = ctbp.tileset_info(input_file, chromsizes)

    assert "max_width" in tileset_info
    assert tileset_info["header"] == header


def test_bedpe_tiles():
    input_file = op.join(testdir, "sample_data", "isidro.bedpe")
    chromsizes_fn = op.join(testdir, "sample_data", "b37.chrom.sizes")

    chromsizes = cs.chromsizes_as_series(chromsizes_fn)
    tileset_info = ctbp.tileset_info(input_file, chromsizes)

    tiles = ctbp.tiles(input_file, ["x.0.0.0"], chromsizes)
    assert len(tiles) > 0
