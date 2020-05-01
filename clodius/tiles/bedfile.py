import math
import pandas as pd
import random

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

    return {
        "max_width": max_width,
        "max_zoom": int(math.log(max_width) / math.log(2)),
        "chromsizes": chromsizes_list,
        "min_pos": [0],
        "max_pos": [max_width]
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

def single_tile(filename, chromsizes, tsinfo, z, x):
    t = pd.read_csv(filename, header=None, delimiter='\t')
    orig_columns = t.columns

    css = chromsizes.cumsum().shift().fillna(0).to_dict()

    t[3] = t[0].map(lambda x: css[x])
    t['xStart'] = t[3] + t[1]
    t['xEnd'] = t[3] + t[2]
    t['ix'] = t.index

    tileStart = x * tsinfo['max_width'] / 2**z
    tileEnd = (x+1) * tsinfo['max_width'] / 2**z

    t = t[t['xEnd'] >= tileStart]
    t = t[t['xStart'] <= tileEnd]

    MAX_PER_TILE = 128

    t = t.sample(MAX_PER_TILE) if len(t) > MAX_PER_TILE else t

    def row_to_bedlike(row):
        ret = {
            'uid': row['ix'],
            'xStart': row['xStart'],
            'xEnd': row['xEnd'],
            'chrOffset': css[row[0]],
            'importance': random.random(),
            'fields': [r for r in row[orig_columns]]
        }

        return ret

    ret = t.apply(row_to_bedlike, axis=1)
    return list(ret.values)

def tiles(filename, tile_ids, chromsizes_map, chromsizes=None):
    tsinfo = tileset_info(filename, chromsizes)

    tile_values = []

    for tile_id in tile_ids:
        tile_option_parts = tile_id.split("|")[1:]
        tile_no_options = tile_id.split("|")[0]
        tile_id_parts = tile_no_options.split(".")
        tile_position = list(map(int, tile_id_parts[1:3]))
        tile_options = dict([o.split(":") for o in tile_option_parts])

        if len(tile_position) < 2:
            raise IndexError("Not enough tile info present")

        chromsizes_id = None
        if "cos" in tile_options:
            chromsizes_id = tile_options["cos"]
        if chromsizes_id in chromsizes_map:
            chromsizes_to_use = chromsizes_map[chromsizes_id]
        else:
            chromsizes_to_use = None

        if not chromsizes_to_use:
            chromsizes_to_use = chromsizes

        z = tile_position[0]
        x = tile_position[1]

        values = single_tile(filename, chromsizes, tsinfo, z, x)

        tile_values += [(tile_id, values)]

    return tile_values
