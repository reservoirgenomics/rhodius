import os.path as op

import clodius.chromosomes as cc
import clodius.tiles.bedfile as ctb


def test_bed_tiles():
	valid_filename = op.join('data', 'regions.valid.bed')
	invalid_filename = op.join('data', 'regions.spaces.bed')

	chromsizes_fn = op.join('data', 'chm13v1.chrom.sizes')

	chromsizes = cc.chromsizes_as_series(chromsizes_fn)
	tiles = ctb.tiles(valid_filename, ['x.0.0'], chromsizes, index_filename=None)

	assert len(tiles) > 0

	tiles = ctb.tiles(invalid_filename, ['x.0.0'], chromsizes, index_filename=None)

	assert 'error' in tiles[0][1]

def test_bed_regions():
	valid_filename = op.join('data', 'regions.valid.bed')
	chromsizes_fn = op.join('data', 'chm13v1.chrom.sizes')
	chromsizes = cc.chromsizes_as_series(chromsizes_fn)

	regions = ctb.regions(valid_filename, chromsizes, 0, 10)
	assert len(regions[0]) == 10

def test_no_item_rgb():
	chromsizes_fn = op.join('data', 'chm13v1.chrom.sizes')
	chromsizes = cc.chromsizes_as_series(chromsizes_fn)
	filename = op.join('data', 'no_item_rgb.bed')

	tiles = ctb.tiles(filename, ['x.0.0'], chromsizes, index_filename=None)
