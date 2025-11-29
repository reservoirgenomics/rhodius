import polars as pl
import pytest
from clodius.models.gff_models import *
from clodius.tiles.gff import parse_gff_to_models

def test_load_and_parse_gff_positions():
    """Test loading positions 1002988 to 1053215 for contig NC_004354.4 from genomic_10k.gff"""
    
    # Load the GFF file
    gff_file = "data/genomic.10k.gff"
    
    # Read GFF file, filtering for the specified contig and position range
    df = pl.read_csv(
        gff_file,
        separator='\t',
        comment_char='#',
        has_header=False,
        new_columns=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    )
    
    # Filter for NC_004354.4 contig and position range 1002988-1053215
    filtered_df = df.filter(
        (pl.col('seqid') == 'NC_004354.4') & 
        (pl.col('start') >= 1002988) & 
        (pl.col('end') <= 1053215)
    )
    
    assert len(filtered_df) > 0, "No entries found in the specified range"
    
    # Parse the filtered dataframe into models
    genes, transcripts = parse_gff_to_models(filtered_df)
    
    # Assertions to verify parsing worked correctly
    assert len(genes) > 0, "No genes were parsed"
    assert len(transcripts) > 0, "No transcripts were parsed"
    
    # Check specific gene exists (TfIIA-S-2 gene should be in this range)
    tfiia_gene = None
    for gene_model in genes.values():
        if 'TfIIA-S-2' in gene_model.gene.attributes.get('Name', ''):
            tfiia_gene = gene_model
            break
    
    assert tfiia_gene is not None, "TfIIA-S-2 gene not found in parsed results"
    assert tfiia_gene.gene.start == 1011101, "TfIIA-S-2 gene start position incorrect"
    assert tfiia_gene.gene.end == 1011643, "TfIIA-S-2 gene end position incorrect"
    assert len(tfiia_gene.transcripts) > 0, "TfIIA-S-2 gene should have transcripts"
    
    # Check that transcripts have exons
    transcript_with_exons = None
    for transcript in transcripts.values():
        if len(transcript.exons) > 0:
            transcript_with_exons = transcript
            break
    
    assert transcript_with_exons is not None, "No transcripts with exons found"
    assert len(transcript_with_exons.exons) >= 1, "Transcript should have at least one exon"
    
    # Check that mRNA transcripts have CDS
    mrna_with_cds = None
    for transcript in transcripts.values():
        if isinstance(transcript, mRNA) and len(transcript.cds) > 0:
            mrna_with_cds = transcript
            break
    
    assert mrna_with_cds is not None, "No mRNA transcripts with CDS found"
    assert len(mrna_with_cds.cds) >= 1, "mRNA transcript should have at least one CDS"
    
    return genes, transcripts