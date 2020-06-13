import clodius.chromosomes as cs
import clodius.tiles.bedpe as ctbp
import os.path as op

testdir = op.realpath(op.dirname(__file__))

import numpy as np


def test_bedpe_tileset_info():
    input_file = op.join(testdir, "sample_data", "isidro.bedpe")
    chromsizes_fn = op.join(testdir, "sample_data", "b37.chrom.sizes")

    chromsizes = cs.chromsizes_as_series(chromsizes_fn)
    tileset_info = ctbp.tileset_info(input_file, chromsizes)

    print("tileset_info:", tileset_info)
