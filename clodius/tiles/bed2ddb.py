import collections as col
import sqlite3

from .utils import tiles_wrapper_2d


def tileset_info(filepath):
    conn = sqlite3.connect(filepath)
    c = conn.cursor()

    row = c.execute("SELECT * from tileset_info").fetchone()
    tileset_info = {
        "zoom_step": row[0],
        "max_length": row[1],
        "assembly": row[2],
        "chrom_names": row[3],
        "chrom_sizes": row[4],
        "tile_size": row[5],
        "max_zoom": row[6],
        "max_width": row[7],
        "min_pos": [1, 1],
        "max_pos": [row[1], row[1]],
    }
    conn.close()

    return tileset_info


# Deprecated. Use `tileset_info()`
def get_2d_tileset_info(filepath):
    return tileset_info(filepath)


def tiles(filepath, tile_ids):
    if len(tile_ids) == 0:
        return []

    is_1d = len(list(tile_ids)[0].split(".")) < 4

    if is_1d:
        return tiles_1d(filepath, tile_ids)

    return tiles_2d(filepath, tile_ids)


def tiles_1d(filepath, tile_ids):
    """
    Generate 1D tiles from this dataset.
    Parameters
    ----------
    filepath: str
        The filename of the sqlite db file
    tile_ids: [str...]
        A list of tile ids of the form
    Returns
    -------
    tiles: [(tile_id, tile_value),...]
        A list of values indexed by the tile position
    """
    to_return = []

    for tile_id in tile_ids:
        _, z, x = tile_id.split(".")
        to_return += [(tile_id, get_1d_tiles(filepath, int(z), int(x))[int(x)])]

    return to_return


def get_1d_tiles(filepath, zoom: int, tile_x_pos: int, num_tiles: int = 1):
    """
    Retrieve a contiguous set of tiles from a 2D db tile file.

    Parameters
    ----------
    filepath: str
        The filename of the sqlite db file
    zoom: int
        The zoom level
    tile_x_pos: int
        The x position of the first tile
    num_tiles: int
        The number of tiles to retrieve

    Returns
    -------
    tiles: {pos: tile_value}
        A set of tiles, indexed by position
    """
    ts_info = tileset_info(filepath)

    conn = sqlite3.connect(filepath)
    c = conn.cursor()

    tile_width = ts_info["max_width"] / 2 ** zoom

    tile_x_start_pos = tile_width * tile_x_pos
    tile_x_end_pos = tile_x_start_pos + (tile_width * num_tiles)

    query = f"""
    SELECT fromX, toX, fromY, toY, chrOffset, importance, fields, uid
    FROM intervals, position_index
    WHERE
        intervals.id=position_index.id AND
        zoomLevel <= {zoom} AND
        rToX >= {tile_x_start_pos} AND
        rFromX <= {tile_x_end_pos}
    UNION
    SELECT fromX, toX, fromY, toY, chrOffset, importance, fields, uid
    FROM intervals, position_index
    WHERE
        intervals.id=position_index.id AND
        zoomLevel <= {zoom} AND
        rToY >= {tile_x_start_pos} AND
        rFromY <= {tile_x_end_pos}
    """

    rows = c.execute(query).fetchall()

    new_rows = col.defaultdict(list)

    for r in rows:
        try:
            uid = r[7].decode("utf-8")
        except AttributeError:
            uid = r[7]

        x_start = r[0]
        x_end = r[1]
        y_start = r[2]
        y_end = r[3]

        for i in range(tile_x_pos, tile_x_pos + num_tiles):
            tile_x_start = i * tile_width
            tile_x_end = (i + 1) * tile_width

            if x_start < tile_x_end and x_end >= tile_x_start:
                # add the position offset to the returned values
                new_rows[i] += [
                    {
                        "xStart": x_start,
                        "xEnd": x_end,
                        "yStart": y_start,
                        "yEnd": y_end,
                        "chrOffset": r[4],
                        "importance": r[5],
                        "uid": uid,
                        "fields": r[6].split("\t"),
                    }
                ]
    conn.close()

    return new_rows


def get_1D_tiles(*args):
    return get_1d_tiles(*args)


def tiles_2d(filepath, tile_ids):
    return tiles_wrapper_2d(
        tile_ids, lambda z, x, y: get_2d_tiles(filepath, z, x, y)[(x, y)]
    )


def get_2d_tiles(db_file, zoom, tile_x_pos, tile_y_pos, numx=1, numy=1):
    """
    Retrieve a contiguous set of tiles from a 2D db tile file.

    Parameters
    ----------
    db_file: str
        The filename of the sqlite db file
    zoom: int
        The zoom level
    tile_x_pos: int
        The x position of the first tile
    tile_y_pos: int
        The y position of the first tile
    numx: int
        The width of the block of tiles to retrieve
    numy: int
        The height of the block of tiles to retrieve

    Returns
    -------
    tiles: {pos: tile_value}
        A set of tiles, indexed by position
    """
    tileset_info = get_2d_tileset_info(db_file)

    conn = sqlite3.connect(db_file)

    c = conn.cursor()
    tile_width = tileset_info["max_width"] / 2 ** zoom

    tile_x_start_pos = tile_width * tile_x_pos
    tile_x_end_pos = tile_x_start_pos + (numx * tile_width)

    tile_y_start_pos = tile_width * tile_y_pos
    tile_y_end_pos = tile_y_start_pos + (numy * tile_width)

    query = """
    SELECT fromX, toX, fromY, toY, chrOffset, importance, fields, uid
    FROM intervals,position_index
    WHERE
        intervals.id=position_index.id AND
        zoomLevel <= {} AND
        rToX >= {} AND
        rFromX <= {} AND
        rToY >= {} AND
        rFromY <= {}
    """.format(
        zoom, tile_x_start_pos, tile_x_end_pos, tile_y_start_pos, tile_y_end_pos
    )

    rows = c.execute(query).fetchall()

    new_rows = col.defaultdict(list)

    for r in rows:
        try:
            uid = r[7].decode("utf-8")
        except AttributeError:
            uid = r[7]

        x_start = r[0]
        x_end = r[1]
        y_start = r[2]
        y_end = r[3]

        for i in range(tile_x_pos, tile_x_pos + numx):
            for j in range(tile_y_pos, tile_y_pos + numy):
                tile_x_start = i * tile_width
                tile_x_end = (i + 1) * tile_width

                tile_y_start = j * tile_width
                tile_y_end = (j + 1) * tile_width

                if (
                    x_start < tile_x_end
                    and x_end >= tile_x_start
                    and y_start < tile_y_end
                    and y_end >= tile_y_start
                ):
                    # add the position offset to the returned values
                    new_rows[(i, j)] += [
                        {
                            "xStart": r[0],
                            "xEnd": r[1],
                            "yStart": r[2],
                            "yEnd": r[3],
                            "chrOffset": r[4],
                            "importance": r[5],
                            "uid": uid,
                            "fields": r[6].split("\t"),
                        }
                    ]
    conn.close()

    return new_rows


def get_2D_tiles(*args):
    return get_2d_tiles(*args)
