import hashlib
import random
import h5py

from clodius.tiles.utils import (
    calc_max_width,
    interval_to_chrom_tiles,
    genome_tile_to_intervals,
)

import math
from clodius.tiles.utils import TilesetInfo


def tileset_info(filename, chromsizes):
    max_zoom = math.ceil(math.log(sum(chromsizes)) / math.log(2))
    max_width = 2**max_zoom

    chromsizes_list = [[chrom, int(size)] for chrom, size in chromsizes.iteritems()]

    with h5py.File(filename, "r") as f:
        max_per_tile = f["info"].attrs["max_per_tile"]

    return {
        "max_width": max_width,
        "max_zoom": int(max_zoom),
        "chromsizes": chromsizes_list,
        "min_pos": [0],
        "max_pos": [max_width],
        "max_per_tile": int(max_per_tile),
    }


def single_chromosome_tile(
    filename, chromsizes, tsinfo: dict, chrom: str, z: int, x: int
):
    f = h5py.File(filename, "r")
    max_per_tile = tsinfo["max_per_tile"]
    css = chromsizes.cumsum().shift().fillna(0).to_dict()
    chrom_len = chromsizes.to_dict()[chrom]

    #     print("x", x)
    max_width = calc_max_width(chrom_len)
    tile_width = max_width / 2**z
    tile_start = tile_width * x
    tile_end = tile_width * (x + 1)
    #     print('chromsizes:', chromsizes)
    #     print("max_per_tile:", max_per_tile)
    #     print("sct", z, x)

    items = []
    tile_pos = x

    if str(z) in f:
        # If the requested zoom is higher than the max then we just
        # return the next lowest zoom
        items += [
            (z, x)
            for x in list(
                f[str(z)][tile_pos * max_per_tile : (tile_pos + 1) * max_per_tile]
            )
            if len(x)
        ]

    while z > 0:
        z -= 1
        tile_pos //= 2

        if str(z) in f:
            # If the requested zoom is higher than the max then we just
            # return the next lowest zoom
            items += [
                (z, x)
                for x in list(
                    f[str(z)][tile_pos * max_per_tile : (tile_pos + 1) * max_per_tile]
                )
                if len(x)
            ]
    #             print("z", z, "tile_pos", tile_pos, items)
    formatted = []

    for z, row in items:
        row = row.decode("utf8")
        parts = row.split("\t")

        start = int(parts[1])
        end = int(parts[2])

        if not (end > tile_start and start < tile_end):
            #             print("ts", tile_start, tile_end)
            #             print("se", start, end)
            #             print("no intersection")
            # doesn't intersect tile
            continue

        ret = {
            "uid": hashlib.md5(row.encode("utf-8")).hexdigest(),
            "zoom": z,
            "xStart": css[parts[0]] + int(parts[1]),
            "xEnd": css[parts[0]] + int(parts[2]),
            "chrOffset": css[parts[0]],
            "importance": random.random(),
            "fields": parts,
        }
        formatted += [ret]

    sorted_formatted = sorted(formatted, key=lambda x: (x["zoom"], x["uid"]))
    return sorted_formatted[:max_per_tile]


def single_genome_tile(filename, chromsizes, tsinfo, z, x):
    chrom_tile_poss = []
    intervals = genome_tile_to_intervals(
        filename, chromsizes, TilesetInfo.parse_obj(tsinfo), z, x
    )
    for interval in intervals:
        if interval[0] >= len(tsinfo["chromsizes"]):
            break
        chrom_name = tsinfo["chromsizes"][interval[0]][0]
        chrom_size = tsinfo["chromsizes"][interval[0]][1]

        chrom_tile_poss += [
            (chrom_name, cz, cx)
            for (cz, cx) in interval_to_chrom_tiles(
                interval[1], interval[2], chrom_size
            )
        ]

    #     print("chrom_tile_poss", chrom_tile_poss)
    chrom_tiles = []
    for chrom, cz, cx in chrom_tile_poss:
        if chrom == "chr10":
            chrom_tiles += single_chromosome_tile(
                filename, chromsizes, tsinfo, chrom, cz, cx
            )

    return chrom_tiles


# tileset_info(filename, chromsizes)
# tile0 = single_genome_tile(filename, chromsizes, tsinfo, 1, 0)
