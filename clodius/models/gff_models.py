from typing import List, Optional, Union, Literal
from pydantic import BaseModel, Field


class BaseGFFEntity(BaseModel):
    """Base class for all GFF entities"""
    id: str
    start: int
    end: int
    strand: Optional[Literal['+', '-', '.']] = None
    score: Optional[float] = None
    phase: Optional[int] = None
    attributes: Optional[dict] = None


class Exon(BaseGFFEntity):
    """Exon entity - can be child of any transcript type"""
    pass


class CDS(BaseGFFEntity):
    """Coding sequence entity - child of mRNA"""
    pass


class Gene(BaseGFFEntity):
    """Root gene entity"""
    gene_biotype: Optional[str] = None
    pseudo: Optional[bool] = False


class mRNA(BaseGFFEntity):
    """Protein-coding transcript"""
    parent_gene_id: str
    exons: List[Exon] = Field(default_factory=list)
    cds: List[CDS] = Field(default_factory=list)


class lnc_RNA(BaseGFFEntity):
    """Long non-coding RNA transcript"""
    parent_gene_id: str
    exons: List[Exon] = Field(default_factory=list)


class miRNA(BaseGFFEntity):
    """Mature microRNA"""
    parent_transcript_id: str
    exons: List[Exon] = Field(default_factory=list)


class primary_transcript(BaseGFFEntity):
    """Precursor RNA transcript"""
    parent_gene_id: str
    exons: List[Exon] = Field(default_factory=list)
    mirnas: List[miRNA] = Field(default_factory=list)


class antisense_RNA(BaseGFFEntity):
    """Antisense RNA transcript"""
    parent_gene_id: str
    exons: List[Exon] = Field(default_factory=list)


class snoRNA(BaseGFFEntity):
    """Small nucleolar RNA transcript"""
    parent_gene_id: str
    exons: List[Exon] = Field(default_factory=list)


class tRNA(BaseGFFEntity):
    """Transfer RNA transcript"""
    parent_gene_id: str
    exons: List[Exon] = Field(default_factory=list)


class rRNA(BaseGFFEntity):
    """Ribosomal RNA transcript"""
    parent_gene_id: str
    exons: List[Exon] = Field(default_factory=list)


class snRNA(BaseGFFEntity):
    """Small nuclear RNA transcript"""
    parent_gene_id: str
    exons: List[Exon] = Field(default_factory=list)


class SRP_RNA(BaseGFFEntity):
    """Signal recognition particle RNA"""
    parent_gene_id: str
    exons: List[Exon] = Field(default_factory=list)


class RNase_P_RNA(BaseGFFEntity):
    """RNase P RNA transcript"""
    parent_gene_id: str
    exons: List[Exon] = Field(default_factory=list)


class RNase_MRP_RNA(BaseGFFEntity):
    """RNase MRP RNA transcript"""
    parent_gene_id: str
    exons: List[Exon] = Field(default_factory=list)


class Pseudogene(BaseGFFEntity):
    """Non-functional gene copy"""
    pseudo: bool = True
    exons: List[Exon] = Field(default_factory=list)


# Union type for all transcript types
TranscriptType = Union[
    mRNA, lnc_RNA, primary_transcript, antisense_RNA, snoRNA, 
    tRNA, rRNA, snRNA, SRP_RNA, RNase_P_RNA, RNase_MRP_RNA
]


class GeneModel(BaseModel):
    """Complete gene model with all associated transcripts"""
    gene: Gene
    transcripts: List[TranscriptType] = Field(default_factory=list)
    
    class Config:
        arbitrary_types_allowed = True


class PseudogeneModel(BaseModel):
    """Pseudogene model"""
    pseudogene: Pseudogene
    
    class Config:
        arbitrary_types_allowed = True


# Union type for all gene forms
GeneForm = Union[GeneModel, PseudogeneModel]