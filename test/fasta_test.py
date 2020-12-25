import os.path as op

import clodius.tiles.fasta as ctf

fasta_filename = op.join("data", "GCA_000350705.1_Esch_coli_KTE11_V1_genomic.short.fna")
fai_filename = op.join(
    "data", "GCA_000350705.1_Esch_coli_KTE11_V1_genomic.short.fna.fai"
)


def test_tileset_info():
    tsinfo = ctf.tileset_info(fai_filename)

    assert "max_zoom" in tsinfo
    assert "max_width" in tsinfo


def test_tiles():
    tiles = ctf.tiles(fasta_filename, fai_filename, ["x.0.0"])

    tsinfo = ctf.tileset_info(fai_filename)

    # make sure we get as much text as we're expecting
    assert len(tiles[0][1]['seq']) == tsinfo['max_pos'][0]
