# @lru_cache
def csv_sequence_tileset_functions(
    filename,
    tile_functions,
    colname=None,
    colnum=None,
    header=True,
    sep=",",
    refrow=None,
):
    """Read a csv file and return a list of sequences.

    Parameters
    ----------
    filename: string
        The name of the csv file
    tile_functions:
        A function that will take a list of sequences as a parameters
        and return tileset_info and tiles functions
    colname: Optional[str]
        The name of the column containing the sequences.
    colnum: Optional[int]
        The column number of the sequence logo file. 0-based.
        Only used if colname is not provided.
    sep: string
        The separator used in the csv file
    refrow: A row to use as a reference sequence when calculating
        alignments. Should be 1-based
    """
    import pandas as pd

    if not header:
        header = None
    else:
        header = 0

    if not colname and not colnum:
        raise ValueError("No colname or colnum specified")

    df = pd.read_csv(filename, header=header, sep=sep)

    if not colname:
        colname = df.columns[colnum]

    sequences = df[colname].values

    if refrow:
        refseq = sequences[refrow - 1]
    else:
        refseq = None

    tf = tile_functions(sequences, refseq=refseq, values=df.to_dict(orient="records"))

    orig_tsinfo = tf["tileset_info"]()
    # Decorate the tileset info function so that it returns
    # the column names as well.
    tf["tileset_info"] = lambda: {"columns": list(df.columns), **orig_tsinfo}

    return tf
