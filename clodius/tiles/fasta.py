import clodius.tiles.chromsizes as cts
import math

from clodius.tiles.utils import parse_tile_id, TilesetInfo


def tileset_info(fai_filename):
    """The tileset info of a fai file returns the tileset."""
    tsinfo = cts.tileset_info(fai_filename)
    tsinfo["max_zoom"] = math.ceil(math.log(tsinfo["max_pos"][0]) / math.log(2))
    tsinfo["max_width"] = 2 ** tsinfo["max_zoom"]
    return tsinfo


def tiles(fasta_filename, fai_filename, tile_ids):
    tsinfo = tileset_info(fai_filename)
    print("tsinfo", tsinfo)
    tsinfo = TilesetInfo(**tsinfo)
    generated_tiles = []
    print("tsinfo:", tsinfo)

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
