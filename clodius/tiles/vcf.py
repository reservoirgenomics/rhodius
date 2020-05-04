import math
import pandas as pd
import random
from pysam import VariantFile

from clodius.tiles.bigwig import abs2genomic

def tileset_info(filename, chromsizes):
    """

    Return the bounds of this tileset. The bounds should encompass the entire
    width of this dataset.

    So how do we know what those are if we don't know chromsizes? We can assume
    that the file is enormous (e.g. has a width of 4 trillion) and rely on the
    browser to pass in a set of chromsizes
    """

    # do this so that we can serialize the int64s in the numpy array
    chromsizes_list = []

    for chrom, size in chromsizes.iteritems():
        chromsizes_list += [[chrom, int(size)]]

    max_width = sum([c[1] for c in chromsizes_list])
    MAX_TILE_WIDTH = 100000

    return {
        "max_width": max_width,
        "max_zoom": int(math.log(max_width) / math.log(2)),
        "chromsizes": chromsizes_list,
        "min_pos": [0],
        "max_pos": [max_width],
        "max_tile_width": MAX_TILE_WIDTH,
    }

# def tiles_wrapper(array, tile_ids, not_nan_array=None):
#     tile_values = []

#     for tile_id in tile_ids:
#         parts = tile_id.split(".")

#         if len(parts) < 3:
#             raise IndexError("Not enough tile info present")

#         z = int(parts[1])
#         x = int(parts[2])

#         ret_array = tiles(array, z, x, not_nan_array).reshape((-1))

#         tile_values += [(tile_id, ctf.format_dense_tile(ret_array))]

#     return tile_values

def single_tile(filename, index_filename, chromsizes, tsinfo, z, x, max_tile_width):
    tile_width = tsinfo['max_width'] / 2 ** z

    if tile_width > max_tile_width:
        return {
            "error": "Tile too wide"
        }

    vcf = VariantFile(filename, index_filename=index_filename)  # auto-detect input format

    start_pos = x * tsinfo['max_width'] / 2**z
    end_pos = (x+1) * tsinfo['max_width'] / 2**z

    css = chromsizes.cumsum().shift().fillna(0).to_dict()

    cids_starts_ends = list(abs2genomic(chromsizes, start_pos, end_pos))
    ret_vals = []

    for (cid, start, end) in cids_starts_ends:
        chrom = chromsizes.index[cid]
        ret_vals += [
            {
                'uid': r.id,
                'importance': random.random(),
                'xStart': css[chrom] + r.start,
                'xEnd': css[chrom] + r.stop,
                'chrOffset': css[chrom],
                "fields": [r.chrom, r.start, r.stop, str(r)]
            } for r in
            vcf.fetch(str(chrom), int(start), int(end))
        ]

    return ret_vals

def tiles(filename, tile_ids, index_filename, chromsizes, max_tile_width=None):
    tsinfo = tileset_info(filename, chromsizes)

    tile_values = []

    for tile_id in tile_ids:
        tile_option_parts = tile_id.split("|")[1:]
        tile_no_options = tile_id.split("|")[0]
        tile_id_parts = tile_no_options.split(".")
        tile_position = list(map(int, tile_id_parts[1:3]))

        if len(tile_position) < 2:
            raise IndexError("Not enough tile info present")

        tile_width = tsinfo["max_width"] / 2 ** int(tile_position[0])

        if max_tile_width and tile_width >= max_tile_width:
            # this tile is larger than the max allowed
            return [
                (
                    tile_id,
                    {
                        "error": f"Tile too large, no data returned. Max tile size: {max_tile_width}"
                    },
                )
            ]

        z = tile_position[0]
        x = tile_position[1]

        values = single_tile(
            filename, index_filename, chromsizes, tsinfo, z, x, max_tile_width
        )

        tile_values += [(tile_id, values)]

    return tile_values
