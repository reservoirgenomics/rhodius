import subprocess
from collections import Counter
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Literal, Optional
import tempfile
import os

DNA_ALPHABET = ["A", "C", "G", "T"]
PROTEIN_ALPHABET = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
]


def run_clustal_omega(sequences, seq_ids=None, seqtype="dna"):
    """
    Align sequences with Clustal Omega.

    Args:
        sequences (list of str): Input nucleotide sequences (unaligned).
        seq_ids (list of str, optional): IDs for sequences (default: numbered).

    Returns:
        alignment (MultipleSeqAlignment): Biopython alignment object.
    """
    if seq_ids is None:
        seq_ids = [f"seq{i}" for i in range(len(sequences))]

    # Create temp fasta input file
    with tempfile.NamedTemporaryFile("w", delete=False) as fasta_file:
        input_fasta = fasta_file.name
        records = [
            SeqRecord(Seq(seq), id=sid, description="")
            for seq, sid in zip(sequences, seq_ids)
        ]
        SeqIO.write(records, fasta_file, "fasta")

    # Create temp output file
    output_fasta = input_fasta + "_aligned.fasta"

    params = [
        "clustalo",
        "-i",
        input_fasta,
        "-o",
        output_fasta,
        "--force",
        "--outfmt=fasta",
    ]

    if seqtype is not None:
        params += [f"--seqtype={seqtype.upper()}"]

    # Run Clustal Omega
    subprocess.run(
        params,
        check=True,
    )

    # Parse alignment
    alignment = AlignIO.read(output_fasta, "fasta")

    # Clean up
    os.remove(input_fasta)
    os.remove(output_fasta)

    return alignment


def make_pwm_from_alignment(
    alignment, pseudocount=0, seqtype: Literal["dna", "protein"] = "dna"
):
    """
    Build a PWM from an aligned set of sequences.

    Args:
        alignment (MultipleSeqAlignment): Aligned sequences.
        pseudocount (int): Pseudocount for smoothing.

    Returns:
        pwm (dict): Dictionary of {base: [probabilities]}.
    """
    if seqtype is None:
        seqtype = "dna" if len(alignment[0]) == 4 else "protein"

    alphabet = DNA_ALPHABET if seqtype == "dna" else PROTEIN_ALPHABET

    seq_length = alignment.get_alignment_length()
    pwm = {base: [] for base in alphabet}
    num_seqs = len(alignment)

    for pos in range(seq_length):
        column = [record.seq[pos] for record in alignment]
        counts = Counter(column)

        total = sum(counts.get(base, 0) + pseudocount for base in alphabet)
        for base in alphabet:
            prob = (counts.get(base, 0) + pseudocount) / total
            pwm[base].append(prob)

    # print(pwm.keys())
    # # print(pwm["M"])
    # for key in pwm:
    #     print(f"{key} {pwm[key][0]:.2} {pwm[key][1]:.2}")
    return pwm


def generate_pwm_from_sequences(
    sequences,
    seq_ids=None,
    pseudocount=0,
    seqtype: Optional[Literal["dna", "protein"]] = None,
):
    """
    Align sequences using Clustal Omega and generate a position weight matrix (PWM).

    Args:
        sequences (list of str): Input nucleotide sequences (unaligned).
        seq_ids (list of str, optional): Sequence IDs (default: numbered seq0, seq1, ...).
        pseudocount (int, optional): Pseudocount for PWM smoothing.

    Returns:
        alignment (MultipleSeqAlignment): Aligned sequences.
        pwm (dict): PWM as dictionary {base: [probabilities]}.
    """
    alignment = run_clustal_omega(sequences, seq_ids, seqtype=seqtype)
    pwm = make_pwm_from_alignment(alignment, pseudocount=pseudocount, seqtype=seqtype)
    return pwm
