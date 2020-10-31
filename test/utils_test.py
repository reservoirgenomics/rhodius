from __future__ import print_function

import os.path as op

import clodius.utils as cu


def test_infer_filetype():
	assert cu.infer_filetype('blah.gff') == 'gff'
	assert cu.infer_filetype('blah.gff.gz') == 'gff'
	assert cu.infer_filetype('blah.xyz') == None
	assert cu.infer_filetype('blah.bam') == 'bam'
	assert cu.infer_filetype('blah.bed.bgz') == 'bedfile'
	assert cu.infer_filetype('blah.bed') == 'bedfile'

def test_infer_datatype():
	assert cu.infer_datatype('gff') == 'bedlike'
	assert cu.infer_datatype('cooler') == 'matrix'
	assert cu.infer_datatype('bedfile') == 'bedlike'
	assert cu.infer_datatype('bam') == 'reads'
