import os.path as op

FILETYPES = {
    "bam": {
        "description": "Read mappings",
        "extensions": [".bam"],
        "datatypes": ["reads", "alignments"]
    },
    "cooler": {
        "description": "multi-resolution cooler file",
        "extensions": [".mcool"],
        "datatypes": ["matrix"],
    },
    "bigwig": {
        "description": "Genomics focused multi-resolution vector file",
        "extensions": [".bw", ".bigwig"],
        "datatypes": ["vector"],
    },
    "bedfile": {
        "description": "BED file",
        "extensions": [".bed", ".bed.gz", ".bed.bgz"],
        "datatypes": ["bedlike", "gene-annotations"],
    },
    "beddb": {
        "description": "SQLite-based multi-resolution annotation file",
        "extensions": [".beddb", ".multires.db"],
        "datatypes": ["bedlike", "gene-annotations"],
    },
    "gff": {
        "description": "General feature format",
        "extensions": [".gff", ".gff.gz", ".gff.bgz"],
        "datatypes": ["bedlike"],
    },
    "hitile": {
        "description": "Multi-resolution vector file",
        "extensions": [".hitile"],
        "datatypes": ["vector"],
    },
    "multivec": {
        "description": "Multi-sample vector file",
        "extensions": [".multivec"],
        "datatypes": ["multivec"],
    },
    "time-interval-json": {
        "description": "Time interval notation",
        "extensions": [".htime"],
        "datatypes": ["time-interval"],
    },
}


def infer_filetype(filename):
    for filetype, meta in FILETYPES.items():
        for ext in meta['extensions']:
            if filename.endswith(ext):
                return filetype

    return None


def infer_datatype(filetype):
    if filetype in FILETYPES:
        return FILETYPES[filetype]["datatypes"][0]

    return None
