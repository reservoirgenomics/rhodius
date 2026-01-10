from Bio import Align
from clodius.alignment import align_sequences, alignment_to_subs, order_by_clustering
from clodius.tiles.csv import csv_sequence_tileset_functions


def calc_chr_offset(chromsizes, chrom_id):
    sum = 0
    for chrom in chromsizes:
        if chrom[0] == chrom_id:
            return sum
        sum += chrom[1]


def align_sequences(seq1, seq2):
    """Align two sequences to each other and return an alignment object."""
    aligner = Align.PairwiseAligner()

    aligner.match_score = 1
    aligner.mismatch_score = -4
    aligner.open_gap_score = -6
    aligner.extend_gap_score = -1

    alignments = aligner.align(seq1, seq2)

    best_alignment = alignments[0]

    return best_alignment


def tile_functions(seqs, refseqs, cluster=None, values=None, chromsizes=None):
    """Return a dictionary of tile functions for the pileup track."""
    longest_seq = sum([c[1] for c in chromsizes])

    def tileset_info():
        return {
            "tile_size": longest_seq,
            "resolutions": [1],
            "max_tile_width": longest_seq,
            "format": "subs",
            "min_pos": [0],
            "max_pos": [longest_seq],
            "chromsizes": chromsizes,
        }

    if cluster == "linkage":
        seqs = order_by_clustering(seqs)

    tile = []
    for i, seq in enumerate(seqs):
        for refseq in refseqs:
            a = align_sequences(refseq["seq"], seq)
            start, end, subs = alignment_to_subs(a)

            chr_offset = calc_chr_offset(chromsizes, refseq["id"])

            tv = {
                "id": f"r{i}_{refseq['id']}",
                "from": start + chr_offset,
                "to": end + chr_offset,
                "substitutions": subs,
                "color": 0,
            }

            if values:
                tv["extra"] = values[i]

        tile.append(tv)

    def tiles(tile_ids):
        tiles = []

        for tile_id in tile_ids:
            parts = tile_id.split(".")
            z = int(parts[1])
            x = int(parts[2])

            if z != 0 and x != 0:
                # return an empty tile
                tiles += [(tile_id, [])]
            else:
                # return the entire tile
                tiles += [(tile_id, tile)]

        return tiles

    return {"tileset_info": tileset_info, "tiles": tiles}


def csv_tileset_info(filename, *csv_args, **csv_kwargs):
    """Get tileset info for a sequence logo file file from
    a csv file.

    Parameters
    ----------
    filename: string
        The name of the csv file
    colname: Optional[str]
        The name of the column containing the sequences.
    colnum: Optional[int]
        The column number of the sequence logo file. 0-based.
        Only used if colname is not provided.
    header: bool
        Whether to assume that a header is present in the csv file
    sep: string
        The separator used in the csv file
    refrow: A row to use as a reference sequence when calculating
        alignments. Should be 1-based
    """
    tf = csv_sequence_tileset_functions(
        filename, tile_functions=tile_functions, *csv_args, **csv_kwargs
    )
    return tf["tileset_info"]()


def csv_tiles(filename, tile_ids, *csv_args, **csv_kwargs):
    tf = csv_sequence_tileset_functions(
        filename, tile_functions=tile_functions, *csv_args, **csv_kwargs
    )

    return tf["tiles"](tile_ids)


def get_local_tiles(filename, *csv_args, **csv_kwargs):
    """Get local higlass tiles for the provided file."""
    tsinfo = csv_tileset_info(filename, *csv_args, **csv_kwargs)
    max_resolution = max(tsinfo["resolutions"])

    tile_ids = []

    for i, res in enumerate(sorted(tsinfo["resolutions"], key=lambda x: -x)):
        for j in range(0, max_resolution // res):
            tile_ids += [f"x.{i}.{j}"]

    tiles = dict(csv_tiles(filename, tile_ids, *csv_args, **csv_kwargs))

    return {"tilesetInfo": tsinfo, "tiles": tiles}
