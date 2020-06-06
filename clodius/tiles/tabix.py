import collections as col
import gzip
import math
import numpy as np
import struct

def load_tbi_idx(index_filename):
    """Load a reduced version of a tabix index so that we can
    go through it and get a sense of how much data will be
    retrieved by a query."""
    with gzip.open(index_filename, 'rb') as f:
        b = bytearray(f.read())

        [_,_,_,_, n_ref, format, col_seq, col_beg, col_end, meta, skip, l_nm] = struct.unpack("<4ciiiiiiii", b[:36])
        c = 36

        names = [n.decode('ascii') for n in b[c:c+l_nm].split(b'\0')]
        c += l_nm

        indeces = []

        for i in range(n_ref):
            n_bin = struct.unpack("<i", b[c:c+4])[0]
            c += 4
            bins = col.defaultdict(list)
            for j in range(n_bin):
                [bin_no, n_chunk] = struct.unpack("<Ii", b[c:c+8])
                c += 8

                bytes_to_read = n_chunk * 2 * 8
                unpack_str = f"<{2 * n_chunk}Q"
                bins[bin_no] = struct.unpack(unpack_str, b[c:c+bytes_to_read])
                c += bytes_to_read

            n_intv = struct.unpack("<i", b[c:c+4])[0]
            c += 4 + 8 * n_intv

            indeces += [bins]

            return dict(zip(names, indeces))

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def reg2bins(begin, end, n_lvls=5, min_shift=14):
    """
    generate key of bins which may overlap the given region,
    check out https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3042176/
    and https://samtools.github.io/hts-specs/tabix.pdf
    for more information.
    Parameters
    ----------
    begin: int
        chromosome position begin
    end: int
        chromosome position end
    n_lvls: int, optional
        cluster level, for tabix, set to 5
    min_shift: int, optional
        minimum shift, for tabix, set to 14
    Returns
    -------
    generator
    """
    begin, end = begin, end
    t, s = 0, min_shift + (n_lvls << 1) + n_lvls
    for l in range(n_lvls + 1):
        b, e = t + (begin >> s), t + (end >> s)
        n = e - b + 1
        for k in range(b, e + 1):
            yield k
            n += 1
        t += 1 << ((l << 1) + l)
        s -= 3

def est_query_size(index, name, start, end):
    ix = index[name]
    total_size = 0

    for bin in list(reg2bins(start,end)):
        if ix[bin]:
            for chunk in chunks(ix[bin], 2):
                total_size += (chunk[1] >> 16) - (chunk[0] >> 16)
    #             print(bin, chunk, ix[bin], (chunk[1] >> 16) - (chunk[0] >> 16))
    return total_size

