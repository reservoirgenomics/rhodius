import negspy.coordinates as nc
import numpy as np
import pandas as pd


def chromsizes_array_to_series(chromsizes):
    """
    Convert an array of [[chrname, size]...] values to a series
    indexed by chrname with size values
    """
    chrnames = [c[0] for c in chromsizes]
    chrvalues = [c[1] for c in chromsizes]

    return pd.Series(np.array([int(c) for c in chrvalues]), index=chrnames)


def chromsizes_as_series(chromsizes_filename):
    """Load chromosome sizes as a pandas series."""
    chromsizes = []

    with open(chromsizes_filename) as f:
        for line in f:
            chromsizes += [line.strip().split("\t")]

    return chromsizes_array_to_series(chromsizes)


def load_chromsizes(chromsizes_filename, assembly=None):
    """
    Load a set of chromosomes from a file or using an assembly
    identifier. If using just an assembly identifier the chromsizes
    will be loaded from the negspy repository.

    Parameters:
    -----------
    chromsizes_filename: string
        The file containing the tab-delimited chromosome sizes
    assembly: string
        Assembly name (e.g. 'hg19'). Not necessary if a chromsizes_filename is passed in
    """
    if chromsizes_filename is not None:
        chrom_info = nc.get_chrominfo_from_file(chromsizes_filename)
        chrom_names = chrom_info.chrom_order
        chrom_sizes = [chrom_info.chrom_lengths[c] for c in chrom_info.chrom_order]
    else:
        if assembly is None:
            raise ValueError("No assembly or chromsizes specified")

        chrom_info = nc.get_chrominfo(assembly)
        chrom_names = nc.get_chromorder(assembly)
        chrom_sizes = nc.get_chromsizes(assembly)

    return (chrom_info, chrom_names, chrom_sizes)
