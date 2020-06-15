import functools as ft
import hashlib
import math
import os
import pandas as pd
import random

import pysam
import slugid

import clodius.tiles.tabix as ctt

cache = []


class LRUCache:
    def __init__(self, capacity):
        self.capacity = capacity
        self.tm = 0
        self.cache = {}
        self.lru = {}

    def get(self, key):
        if key in self.cache:
            self.lru[key] = self.tm
            self.tm += 1
            return self.cache[key]
        return None

    def set(self, key, value):
        if len(self.cache) >= self.capacity:
            # find the LRU entry
            old_key = min(self.lru.keys(), key=lambda k: self.lru[k])
            self.cache.pop(old_key)
            self.lru.pop(old_key)
        self.cache[key] = value
        self.lru[key] = self.tm
        self.tm += 1


cache = LRUCache(1)


def tileset_info(filename, chromsizes=None):
    """

    Return the bounds of this tileset. The bounds should encompass the entire
    width of this dataset.

    So how do we know what those are if we don't know chromsizes? We can assume
    that the file is enormous (e.g. has a width of 4 trillion) and rely on the
    browser to pass in a set of chromsizes
    """

    # do this so that we can serialize the int64s in the numpy array
    chromsizes_list = []

    t = pd.read_csv(filename, nrows=2, sep="\t", comment="#", header=None)
    t_head = pd.read_csv(
        filename, nrows=2, sep="\t", comment="#", header=None, skiprows=1
    )

    if (t.dtypes == t_head.dtypes).all():
        has_header = False
        header = ""
    else:
        header = "\t".join(t.head().values[0])
        has_header = True

    if chromsizes is None:
        return {"error": "No chromsizes found. Make sure the assembly: tag is set"}

    for chrom, size in chromsizes.iteritems():
        chromsizes_list += [[chrom, int(size)]]

    max_width = sum([c[1] for c in chromsizes_list])

    if os.stat(filename).st_size > 20e6:
        return {"error": "File too large (>20Mb), please index"}

    tsinfo = {
        "max_width": max_width,
        "max_zoom": int(math.log(max_width) / math.log(2)),
        "chromsizes": chromsizes_list,
        "min_pos": [0, 0],
        "max_pos": [max_width, max_width],
        "header": header,
    }

    return tsinfo


def row_to_bedlike(row, css, orig_columns):
    ret = {
        "uid": row["ix"],
        "xStart": row["xStart"],
        "xEnd": row["xEnd"],
        "yStart": row["yStart"],
        "yEnd": row["yEnd"],
        "xChrOffset": css[str(row[0])],
        "yChrOffset": css[str(row[3])],
        "importance": random.random(),
        "fields": [r for r in row[orig_columns]],
    }

    return ret


def get_compression(filename):
    if filename.endswith(".gz..") or filename.endswith(".gz"):
        return "gzip"
    elif filename.endswith(".bz2..") or filename.endswith(".bz2"):
        return "bz2"
    elif filename.endswith(".zip..") or filename.endswith(".zip"):
        return "zip"
    elif filename.endswith(".xz..") or filename.endswith("xz"):
        return "xz"
    else:
        return None


def ts_hash(filename, chromsizes):
    cs_hash = hashlib.md5(str(chromsizes).encode("utf-8")).hexdigest()
    return f"{filename}.{cs_hash}"


def single_tile(filename, chromsizes, tsinfo, z, x, y):
    hash_ = ts_hash(filename, chromsizes)

    # hash the loaded data table so that we don't have to read the entire thing
    # and calculate cumulative start and end positions
    val = cache.get(hash_)

    if val is None:
        skiprows = 0

        # if this file has a header, skip the first row
        if len(tsinfo["header"]):
            skiprows = 1

        t = pd.read_csv(
            filename,
            header=None,
            comment="#",
            sep="\t",
            compression=get_compression(filename),
            skiprows=skiprows,
        )

        cache.set(hash_, t)

        orig_columns = t.columns
        css = chromsizes.cumsum().shift().fillna(0).to_dict()

        # xStart and xEnd are cumulative start and end positions calculated
        # as if the chromosomes are concatenated from end to end
        t["xChromStart"] = [
            css[str(x)] for x in t[0].values
        ]  # .astype("str").map(lambda x: css[str(x)])
        t["yChromStart"] = [
            css[str(x)] for x in t[3].values
        ]  # .astype("str").map(lambda x: css[str(x)])

        t["xStart"] = t["xChromStart"] + t[1]
        t["xEnd"] = t["xChromStart"] + t[2]

        t["yStart"] = t["yChromStart"] + t[4]
        t["yEnd"] = t["yChromStart"] + t[5]

        t["ix"] = t.index

        val = {"rows": t, "orig_columns": orig_columns, "css": css}
        cache.set(hash_, val)

    t = val["rows"]
    orig_columns = val["orig_columns"]
    css = val["css"]

    xTileStart = x * tsinfo["max_width"] / 2 ** z
    xTileEnd = (x + 1) * tsinfo["max_width"] / 2 ** z

    yTileStart = y * tsinfo["max_width"] / 2 ** z
    yTileEnd = (y + 1) * tsinfo["max_width"] / 2 ** z

    t = t.query(
        f"xEnd >= {xTileStart} & xStart <= {xTileEnd} & "
        + f"yEnd >= {yTileStart} & yStart <= {yTileEnd}"
    )
    MAX_PER_TILE = 512

    t = t.sample(MAX_PER_TILE) if len(t) > MAX_PER_TILE else t

    ret = t.apply(
        ft.partial(row_to_bedlike, css=css, orig_columns=orig_columns), axis=1
    )
    return list(ret.values)


def tiles(filename, tile_ids, chromsizes):
    tsinfo = tileset_info(filename, chromsizes)

    tile_values = []

    for tile_id in tile_ids:
        tile_option_parts = tile_id.split("|")[1:]
        tile_no_options = tile_id.split("|")[0]
        tile_id_parts = tile_no_options.split(".")
        tile_position = list(map(int, tile_id_parts[1:4]))
        tile_options = dict([o.split(":") for o in tile_option_parts])

        if len(tile_position) < 3:
            raise IndexError("Not enough tile info present")

        z = tile_position[0]
        x = tile_position[1]
        y = tile_position[2]

        values = single_tile(filename, chromsizes, tsinfo, z, x, y)

        tile_values += [(tile_id, values)]

    return tile_values
