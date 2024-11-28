import json
import math

import numpy as np

import clodius.tiles.bigwig as ctbw
from clodius.tiles.tabix import est_query_size_ix, load_bai_index
from clodius.tiles.utils import abs2genomic
from clodius.utils import TILE_OPTIONS_CHAR

import oxbow as ox
import polars as pl


def get_cigar_substitutions(pos, query_length, cigartuples):
    subs = []
    curr_pos = 0

    cigartuples = cigartuples
    readstart = pos
    readend = pos + query_length

    for ctuple in cigartuples:
        if ctuple[0] == "X":
            subs.append((readstart + curr_pos, "X", ctuple[1]))
            curr_pos += ctuple[1]
        elif ctuple[0] == "I":
            subs.append((readstart + curr_pos, "I", ctuple[1]))
        elif ctuple[0] == "D":
            subs.append((readstart + curr_pos, "D", ctuple[1]))
            curr_pos += ctuple[1]
        elif ctuple[0] == "N":
            subs.append((readstart + curr_pos, "N", ctuple[1]))
            curr_pos += ctuple[1]
        elif ctuple[0] == "M" or ctuple[0] == "=":
            curr_pos += ctuple[1]

    if len(cigartuples):
        first_ctuple = cigartuples[0]
        last_ctuple = cigartuples[-1]

        if first_ctuple[0] == "S":
            subs.append((readstart - first_ctuple[1], "S", first_ctuple[1]))
        if first_ctuple[0] == "H":
            subs.append((readstart - first_ctuple[1], "H", first_ctuple[1]))

        if last_ctuple[0] == "S":
            subs.append((readend - last_ctuple[1], "S", last_ctuple[1]))
        if last_ctuple[0] == "H":
            subs.append((readend, "H", last_ctuple[1]))

    return subs


def parse_cigar_string(cigar):
    parts = []
    curr = 0
    for c in cigar:
        if c.isnumeric():
            curr = curr * 10 + int(c)
        else:
            parts += [(c, curr)]
            curr = 0
    return parts


def reconstruct_ref(seq, md, cigar):
    """Reconstruct a reference sequence that has the insertions from the query sequence.

    The reason we can't exclude the insertions is that they are encoded for in the CIGAR
    string so we would need to use that to remove them.
    """
    ref = ""
    i_seq = 0
    i_md = 0
    match_count = 0
    deletion = False
    l = 0

    new_seq = ""

    # go through the cigar and remove the ignored bases
    for i_cig in range(len(cigar)):
        if cigar[i_cig].isnumeric():
            # getting the number of bases the upcoming operation applies to
            l = l * 10 + int(cigar[i_cig])
        else:
            op = cigar[i_cig]
            if op == "S" or op == "I":
                i_seq += l
            else:
                new_seq += seq[i_seq : i_seq + l]

            l = 0
            i_cig += 1

    i_seq = 0
    seq = new_seq

    # let's iterate over the entire md string
    for i_md in range(len(md)):
        # if we encounter a numeric value then we keep track of what it is
        if md[i_md].isnumeric():
            match_count = match_count * 10 + int(md[i_md])
            # We're definitely no in a deletion if we're in a numeric number
            deletion = False
        else:
            # Add the matches that we've gone over
            # If we've been going over a deletion or mismatches, then match_count will be 0
            ref += seq[i_seq : i_seq + match_count]
            i_seq += match_count
            match_count = 0

            if md[i_md] == "^":
                # We're starting a deletion sequence
                deletion = True
            else:
                # A letter can indicate that we're either encountering a deletion
                # or a mistmatch

                if deletion:
                    # It's a deletion in the reference
                    ref += md[i_md]
                else:
                    # It's a mismatch, add the MD letter and skip the sequence letter
                    ref += md[i_md]
                    i_seq += 1

    # Add the last match_count stretch
    ref += seq[i_seq : i_seq + match_count]
    return ref


def variants_list(seq, ref, pos, cigar):
    """Get a list of variants that are in seq relative to ref

    Returns:
        A list of (query_pos, ref_pos, query_base) pairs.
    """
    l = 0
    i_cig = 0
    i_seq = 0
    i_ref = 0
    variants = []

    for i_cig in range(len(cigar)):
        if cigar[i_cig].isnumeric():
            # getting the number of bases the upcoming operation applies to
            l = l * 10 + int(cigar[i_cig])
        else:
            op = cigar[i_cig]
            # print("l", l, "op", op)
            if op == "M" or op == "=" or op == "X":
                # Match or mismatch
                # print("l:", l, "op", op)
                for j in range(l):
                    # print(
                    #     "comparing",
                    #     i_seq + j,
                    #     seq[i_seq + j],
                    #     i_ref + j,
                    #     ref[i_ref + j],
                    # )

                    if seq[i_seq + j] != ref[i_ref + j]:
                        # print(
                        #     "mm",
                        #     i_seq + j,
                        #     i_ref + j,
                        #     seq[i_seq + j],
                        #     ref_with_insertions[i_ref + j],
                        # )
                        variants += [(i_seq + j, pos + i_ref + j - 1, seq[i_seq + j])]

                i_seq += l
                i_ref += l
            elif op == "I" or op == "S":
                i_seq += l
            elif op == "D":
                # The query has had deletions from the reference added to it so we have to skip
                # them in the reference
                i_ref += l
            l = 0

        i_cig += 1

    return variants


def load_reads(file, start_pos, end_pos, chromsizes=None, index_file=None, cache=None):
    """
    Sample reads from the specified region, assuming that the chromosomes
    are ordered in some fashion. Returns an list of reads

    Parameters:
    -----------
    file: file-like
        The opened BAM file
    start_pos: int
        The start position of the sampled region
    end_pos: int
        The end position of the sampled region
    chromsize: pandas.Series
        A listing of chromosome sizes. If not provided, the chromosome
        list will be extracted from the the bam file header
    index_file: file-like
        The index file
    cache:
        An object that implements the `get`, `set` and `exists` methods
        for caching data

    Returns
    -------
    reads: [read1, read2...]
        The list of in the sampled regions
    """
    # if chromorder is not None...
    # specify the chromosome order for the fetched reads
    if isinstance(file, str):
        file = open(file, "rb")
    if index_file and isinstance(index_file, str):
        index_file = open(index_file, "rb")

    if chromsizes is not None:
        chromsizes_list = []

        for chrom, size in chromsizes.items():
            chromsizes_list += [[chrom, int(size)]]
    else:
        file.seek(0)
        ipc = ox.read_bam_references(file)
        chroms = pl.read_ipc(ipc)
        references = np.array(chroms["name"])
        lengths = np.array(chroms["length"])

        ref_lengths = dict(zip(references, lengths))

        # we're going to create a natural ordering for references
        # e.g. (chr1, chr2,..., chr10, chr11...chr22,chrX, chrY, chrM...)
        references = ctbw.natsorted(references)
        lengths = [ref_lengths[r] for r in references]
        chromsizes_list = list(zip(references, [int(l) for l in lengths]))

    lengths = [r[1] for r in chromsizes_list]

    abs_chrom_offsets = np.r_[0, np.cumsum(lengths)]

    results = {
        "id": [],
        "from": [],
        "to": [],
        "md": [],
        "chrName": [],
        "chrOffset": [],
        "cigar": [],
        "m1From": [],
        "m1To": [],
        "m2From": [],
        "m2To": [],
        "mapq": [],
        "tags.HP": [],
        "strand": [],
        "variants": [],
        "cigars": [],
    }

    strands = {True: "-", False: "+"}
    import time

    index_file.seek(0)
    idx = load_bai_index(index_file)

    total_size = 0
    # check the size of the file to load to get an approximation
    # of whether we're going to return too much data
    for cid, start, end in abs2genomic(lengths, start_pos, end_pos):
        if cid >= len(chromsizes_list):
            continue
        total_size += est_query_size_ix(idx[cid], start, end)

    MAX_SIZE = 4e6
    if total_size > MAX_SIZE:
        return {"error": "Tile encompasses too much data: {total_size}"}

    for cid, start, end in abs2genomic(lengths, start_pos, end_pos):
        chr_offset = int(abs_chrom_offsets[cid])

        if cid >= len(chromsizes_list):
            continue

        seq_name = f"{chromsizes_list[cid][0]}"
        if start == 0:
            start = 1

        file.seek(0)
        index_file.seek(0)
        # ipc = ox.read_bam(
        #     file, f"{seq_name}:{start}-{end}", index=index_file, tags=set(["MD", "HP"])
        # )
        region = f"{seq_name}:{start}-{end}"
        ipc = ox.read_bam(file, region, index=index_file)
        reads_df = pl.read_ipc(ipc)
        # We can drastically speed these functions up by coding them in Rust in oxbow
        results["cigars"] = [
            get_cigar_substitutions(pos - 1, end - pos, parse_cigar_string(cigar))
            for pos, end, cigar in zip(
                reads_df["pos"], reads_df["end"], reads_df["cigar"]
            )
        ]

        num_reads = len(reads_df)

        results["first_seq"] = list(reads_df["flag"] & 64)
        results["last_seq"] = list(reads_df["flag"] & 128)
        results["is_paired"] = list(reads_df["flag"] & 1)
        results["from"] = list(reads_df["pos"] - 1)
        results["to"] = list(reads_df["end"])
        results["chrName"] = list(reads_df["rname"])
        results["chrOffset"] = [chr_offset] * num_reads

        results["id"] = [
            name if not is_paired else (f"{name}_1" if first else f"{name}_2")
            for name, first, is_paired in zip(
                reads_df["qname"], results["first_seq"], results["is_paired"]
            )
        ]
        results["readName"] = list(reads_df["qname"])

        if "HP" not in reads_df:
            results["tags.HP"] = [0] * num_reads
        else:
            results["tags.HP"] = reads_df["HP"]

        if "MD" not in reads_df:
            results["md"] = [""] * num_reads
            results["variants"] = []
        else:
            results["md"] = list(reads_df["MD"])
            results["variants"] = [
                (
                    variants_list(
                        iseq, reconstruct_ref(iseq, imd, icigar), ipos, icigar
                    )
                    if imd
                    else []
                )
                for iseq, imd, ipos, icigar in zip(
                    reads_df["seq"], reads_df["MD"], reads_df["pos"], reads_df["cigar"]
                )
            ]

    return results


def get_cached_variants(cache, read_id):
    """Try to get variants from a read we've seen before.

    This is useful for ONT reads where there's many variants
    per read and retrieving them takes a while.
    """
    cache_id = f"variants.{read_id}"
    if cache and cache.exists(cache_id):
        return json.loads(cache.get(cache_id))

    return None


def set_cached_variants(cache, read_id, variants):
    """Save a set of variants to the cache."""
    cache_id = f"variants.{read_id}"
    if cache:
        cache.set(cache_id, json.dumps(variants))


def alignment_tileset_info(file, chromsizes):
    """
    Get the tileset info for a bam file

    Parameters
    ----------
    tileset: tilesets.models.Tileset object
        The tileset that the tile ids should be retrieved from

    Returns
    -------
    tileset_info: {'min_pos': [],
                    'max_pos': [],
                    'tile_size': 1024,
                    'max_zoom': 7
                    }
    """
    if isinstance(file, str):
        file = open(file, "rb")

    if chromsizes is not None:
        chromsizes_list = []

        for chrom, size in chromsizes.items():
            chromsizes_list += [[chrom, int(size)]]

        total_length = sum([c[1] for c in chromsizes_list])
    else:
        file.seek(0)
        ipc = ox.read_bam_references(file)
        chroms = pl.read_ipc(ipc)
        total_length = sum(chroms["length"])

        references = np.array(chroms["name"])
        lengths = np.array(chroms["length"])

        ref_lengths = dict(zip(references, lengths))
        references = ctbw.natsorted(references)

        lengths = [ref_lengths[r] for r in references]
        chromsizes_list = list(zip(references, [int(l) for l in lengths]))

    tile_size = 256
    max_zoom = math.ceil(math.log(total_length / tile_size) / math.log(2))

    # this should eventually be a configurable option
    MAX_TILE_WIDTH = 100000

    tileset_info = {
        "min_pos": [0],
        "max_pos": [total_length],
        "max_width": tile_size * 2**max_zoom,
        "tile_size": tile_size,
        "chromsizes": chromsizes_list,
        "max_zoom": max_zoom,
        "max_tile_width": MAX_TILE_WIDTH,
    }

    return tileset_info


def alignment_tiles(
    file,
    tile_ids,
    index_file=None,
    chromsizes=None,
    max_tile_width=None,
    cache=None,
):
    """
    Generate tiles from a bigwig file.

    Parameters
    ----------
    tileset: tilesets.models.Tileset object
        The tileset that the tile ids should be retrieved from
    tile_ids: [str,...]
        A list of tile_ids (e.g. xyx.0.0) identifying the tiles
        to be retrieved
    index_filename: str
        The name of the file containing the index
    max_tile_width: int
        How wide can each tile be before we return no data. This
        can be used to limit the amount of data returned.
    cache:
        An object that implements the `get`, `set` and `exists` methods
        for caching data
    Returns
    -------
    tile_list: [(tile_id, tile_data),...]
        A list of tile_id, tile_data tuples
    """
    if index_file is None:
        if isinstance(file, str):
            index_file = file + ".bai"
        else:
            raise ValueError(
                "A file pointer is provided without an index file. "
                "Please specify an index file"
            )
    generated_tiles = []
    tsinfo = alignment_tileset_info(file, chromsizes)

    for tile_id in tile_ids:
        tile_id_parts = tile_id.split(TILE_OPTIONS_CHAR)[0].split(".")
        tile_position = list(map(int, tile_id_parts[1:3]))

        tile_width = tsinfo["max_width"] // 2 ** int(tile_position[0])

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
        else:
            start_pos = int(tile_position[1]) * tile_width
            end_pos = start_pos + tile_width

            tile_value = load_reads(
                file,
                start_pos=start_pos,
                end_pos=end_pos,
                chromsizes=chromsizes,
                index_file=index_file,
                cache=cache,
            )
            generated_tiles += [(tile_id, tile_value)]

    return generated_tiles


def tileset_info(file, chromsizes=None):
    return alignment_tileset_info(file, chromsizes)


def tiles(
    file,
    tile_ids,
    index_file=None,
    chromsizes=None,
    max_tile_width=None,
    cache=None,
):
    return alignment_tiles(
        file,
        tile_ids,
        index_file=index_file,
        chromsizes=chromsizes,
        max_tile_width=None,
        cache=cache,
    )
