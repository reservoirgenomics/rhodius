import math
from typing import Any, List, Tuple

import numpy as np

import clodius.tiles.chromsizes as cts
from clodius.tiles.format import format_dense_tile
from clodius.tiles.utils import TilesetInfo, abs2genome_fn, parse_tile_id
from pysam import FastaFile

TILE_SIZE = 1024

def convert_bases_to_multivec(seq):
    res = []

    to_append = {
        'a': [1, 0, 0, 0, 0, 0],
        't': [0, 1, 0, 0, 0, 0],
        'g': [0, 0, 1, 0, 0, 0],
        'c': [0, 0, 0, 1, 0, 0],
        'n': [0, 0, 0, 0, 1, 0],
    }


    for c in seq:
        res.append(to_append.get(c.lower(), [0, 0, 0, 0, 0, 1]))
    
    return res;

def tileset_info(fai_filename):
    """The tileset info of a fai file returns the tileset."""
    tsinfo = cts.tileset_info(fai_filename)
    tsinfo["max_zoom"] = math.ceil(
        math.log(tsinfo["max_pos"][0] / TILE_SIZE) / math.log(2)
    )

    tsinfo["max_width"] = TILE_SIZE * 2 ** tsinfo["max_zoom"]
    # tsinfo['bins_per_dimension'] = TILE_SIZE
    tsinfo['tile_size'] = TILE_SIZE
    tsinfo['datatype'] = 'multivec_singleres_sequence'
    return tsinfo

def tiles(fasta_filename: str, tile_ids: List[str],
    index_filename: str, chromsizes_fn: str=None) -> List[Tuple[str, Any]]:
    """Retrieve higlass tiles.

    Arguments:
        fasta_filename: The name of the fasta file to load
        tile_ids: The incoming tile ids (e.g. 'x.0.0')
        fai_filename: The name of the fasta index file (`samtools faidx`)
        chromsizes_filename: The chromsizes filename to use in case we
            want a specific chromosome order.
    Returns:
        Tile data
    """
    tsinfo = tileset_info(index_filename)
    tsinfo = TilesetInfo(**tsinfo)
    generated_tiles = []

    fa_file = FastaFile(fasta_filename, index_filename)

    if not chromsizes_fn:
        chromsizes_fn = index_filename

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
            chromsizes_fn, tile_info.start[0], tile_info.end[0]
        ):
            seq += fa_file.fetch(chr_interval.name, chr_interval.start, chr_interval.end)

        res = convert_bases_to_multivec(seq)
        tile = format_dense_tile(np.array(res).T)
        tile['shape'] = [6, len(seq)]

        generated_tiles += [(tile_id, tile)]

    return generated_tiles
