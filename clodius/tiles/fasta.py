import math

from pydantic import BaseModel

import clodius.tiles.chromsizes as cts
from clodius.tiles.utils import TilesetInfo, abs2genome_fn, parse_tile_id
from pysam import FastaFile

TILE_SIZE = 1024

class FastaTile(BaseModel):
    seq: str

def tileset_info(fai_filename):
    """The tileset info of a fai file returns the tileset."""
    tsinfo = cts.tileset_info(fai_filename)
    tsinfo["max_zoom"] = math.ceil(
        math.log(tsinfo["max_pos"][0] / TILE_SIZE) / math.log(2)
    )
    tsinfo["max_width"] = TILE_SIZE * 2 ** tsinfo["max_zoom"]
    return tsinfo

def tiles(fasta_filename, fai_filename, tile_ids):
    tsinfo = tileset_info(fai_filename)
    tsinfo = TilesetInfo(**tsinfo)
    generated_tiles = []

    fa_file = FastaFile(fasta_filename, fai_filename)

    for tile_id in tile_ids:
        tile_info = parse_tile_id(tile_id, tsinfo)

        zoom_diff = tsinfo.max_zoom - tile_info.zoom
        if zoom_diff > 3:
            generated_tiles += [
                (
                    tile_id,
                    {
                        "error": f"Tile too wide (zoom level {tile_info.zoom}). Please zoom in."
                    },
                )
            ]
            continue

        seq = ""

        for chr_interval in abs2genome_fn(
            fai_filename, tile_info.start[0], tile_info.end[0]
        ):
            seq += fa_file.fetch(chr_interval.name, chr_interval.start, chr_interval.end)

        tile = FastaTile(seq=seq)
        generated_tiles += [(tile_id, tile.dict())]

    return generated_tiles
