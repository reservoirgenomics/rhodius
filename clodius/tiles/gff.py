import functools as ft
import random

import clodius.tiles.bedfile as ctb
import pandas as pd


def gff_chromsizes(filename):
    """Use the "regions" sections of a GFF file as the chromsizes."""
    t = pd.read_csv(
       filename, header=None, delimiter="\t", comment='#', compression=ctb.get_compression(filename)
    )
    regions = t[t[2] == 'region']
    return pd.Series(regions[4].values, index=regions[0])

def row_to_bedlike(row, css, orig_columns):
    attrs = dict([x.split('=') for x in row[8].split(';')])

    ret = {
        "uid": row["ix"],
        "xStart": row["xStart"],
        "xEnd": row["xEnd"],
        "chrOffset": css[row[0]],
        "importance": random.random(),
        "fields": [row[0], row[3], row[4], attrs['Name']],
    }

    return ret

def tileset_info(filename, chromsizes=None, index_filename=None):
    """

    Return the bounds of this tileset. The bounds should encompass the entire
    width of this dataset.

    So how do we know what those are if we don't know chromsizes? We can assume
    that the file is enormous (e.g. has a width of 4 trillion) and rely on the
    browser to pass in a set of chromsizes
    """
    if chromsizes is None:
        chromsizes = gff_chromsizes(filename)

    return ctb.tileset_info(filename, chromsizes, index_filename)

def single_tile(filename, chromsizes, tsinfo, z, x, settings=None):
    hash_ = ctb.ts_hash(filename, chromsizes)

    if settings is None:
        settings = {}
    # hash the loaded data table so that we don't have to read the entire thing
    # and calculate cumulative start and end positions
    val = ctb.cache.get(hash_)
    print("hey")

    if val is None:
        t = pd.read_csv(
            filename, comment='#', header=None, delimiter="\t", compression=ctb.get_compression(filename)
        )
        t = t[t[2] == 'gene']

        orig_columns = t.columns
        css = chromsizes.cumsum().shift().fillna(0).to_dict()

        # xStart and xEnd are cumulative start and end positions calculated
        # as if the chromosomes are concatenated from end to end
        t["chromStart"] = t[0].map(lambda x: css[x])
        t["xStart"] = t["chromStart"] + t[3]
        t["xEnd"] = t["chromStart"] + t[4]
        t["ix"] = t.index

        val = {"rows": t, "orig_columns": orig_columns, "css": css}
        ctb.cache.set(hash_, val)

    t = val["rows"]
    orig_columns = val["orig_columns"]
    css = val["css"]

    tileStart = x * tsinfo["max_width"] / 2 ** z
    tileEnd = (x + 1) * tsinfo["max_width"] / 2 ** z

    t = t.query(f"xEnd >= {tileStart} & xStart <= {tileEnd}")
    MAX_PER_TILE = settings.get("MAX_BEDFILE_ENTRIES") or 1024

    t = t.sample(MAX_PER_TILE) if len(t) > MAX_PER_TILE else t

    ret = t.apply(
        ft.partial(row_to_bedlike, css=css, orig_columns=orig_columns), axis=1
    )
    return list(ret.values)

def tiles(filename, tile_ids, chromsizes=None, index_filename=None, settings=None):
    if chromsizes is None:
        chromsizes = gff_chromsizes(filename)

    return ctb.tiles(filename, tile_ids, chromsizes, index_filename, settings, single_tile_func = single_tile)
