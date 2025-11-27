import functools as ft
import random

import clodius.tiles.bedfile as ctb
import pandas as pd

from clodius.utils import get_file_compression


def gff_chromsizes(filename):
    """Use the "regions" sections of a GFF file as the chromsizes."""
    if isinstance(filename, str):
        filename = open(filename, "rb")

    t = pd.read_csv(
        filename,
        header=None,
        delimiter="\t",
        comment="#",
        compression=get_file_compression(filename),
    )
    regions = t[t[2] == "region"]
    return pd.Series(regions[4].values, index=regions[0])


def row_to_bedlike(row, css, orig_columns):
    attrs = dict([x.split("=") for x in row[8].split(";")])

    ret = {
        "uid": row["ix"],
        "xStart": row["xStart"],
        "xEnd": row["xEnd"],
        "chrOffset": css[row[0]],
        "importance": random.random(),
        "fields": [row[0], row[3], row[4], attrs["Name"], "-", row[6]],
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


def rows_to_genes(rows, css):
    """Convert a set of gff rows into gene annotations in BED12+3 format.

    From https://genome.ucsc.edu/FAQ/FAQformat.html#format1.7

    The format consists of the following:

    chrom - Name of the chromosome (or contig, scaffold, etc.).
    chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
    chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
    name - Name given to a region (preferably unique). Use "." if no name is assigned.
    score - Indicates how dark the peak will be displayed in the browser (0-1000). If all scores were "0" when the data were submitted to the DCC, the DCC assigned scores 1-1000 based on signal value. Ideally the average signalValue per base spread is between 100-1000.
    strand - +/- to denote strand or orientation (whenever applicable). Use "." if no orientation is assigned.
    thickStart - The starting position at which the feature is drawn thickly. Not used in gappedPeak type, set to 0.
    thickEnd - The ending position at which the feature is drawn thickly. Not used in gappedPeak type, set to 0.
    itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). Not used in gappedPeak type, set to 0.
    blockCount - The number of blocks (exons) in the BED line.
    blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
    blockStarts - A comma-separated list of block starts. The first value must be 0 and all of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.
    signalValue - Measurement of overall (usually, average) enrichment for the region.
    pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
    qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
    """
    # GFF file entries may have PARENT=<ID> hierarchies

    pass


def single_tile(filename, chromsizes, tsinfo, z, x, settings=None):
    if isinstance(filename, str):
        filename = open(filename, "rb")

    hash_ = ctb.ts_hash(filename, chromsizes)

    if settings is None:
        settings = {}
    # hash the loaded data table so that we don't have to read the entire thing
    # and calculate cumulative start and end positions
    val = ctb.cache.get(hash_)

    if val is None:
        t = pd.read_csv(
            filename,
            comment="#",
            header=None,
            delimiter="\t",
            compression=get_file_compression(filename),
        )
        t = t[t[2] == "gene"]

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

    tileStart = x * tsinfo["max_width"] / 2**z
    tileEnd = (x + 1) * tsinfo["max_width"] / 2**z

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

    return ctb.tiles(
        filename,
        tile_ids,
        chromsizes,
        index_filename,
        settings,
        single_tile_func=single_tile,
    )
