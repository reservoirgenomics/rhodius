import os.path as op

import clodius.tiles.chromsizes as ctcs


def test_get_tileset_info():
    filename = op.join("data", "hg38.chrom.sizes")

    tsinfo = ctcs.tileset_info(filename)
    assert len(tsinfo['chromsizes']) > 2
    # TODO: Do something with the return value
