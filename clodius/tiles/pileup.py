from Bio import Align
from Bio.Seq import Seq

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
    
def get_subs(a):
    """Convert a BioPython alignment object into the pileup "subs" format.
    
    This format lists the modifications to the target sequence that turn it into
    the query.
    """
    parts = []
    ttrue = 0
    tpos = 0
    qpos = 0

    start = 0
    end = 0

    aligneds = list(zip(a.aligned[0], a.aligned[1]))

    for i, ((ts, te), (qs, qe)) in enumerate(aligneds):
        ts,te,qs,qe = int(ts), int(te), int(qs), int(qe)
        
        if i == 0:
            # start position
            start = ts
            tpos = ts
            ttrue = 0
        if i == len(aligneds) - 1:
            # end position
            end = te
            
        if ts > tpos:
            parts += [{'pos': ttrue, 'type': 'D', 'length': ts - tpos}]
            ttrue += ts - tpos
        if qs > qpos:
            parts += [{'pos': ttrue, 'type': 'I', 'length': qs - qpos}]
        for i in range(te - ts):
            if a.target[ts + i] != a.query[qs + i]:
                parts += [{'pos': ttrue + i, 'type': 'X', 'length': 1, 'base': a.target[ts + i], 'variant': a.query[qs + i]}]

        ttrue += (te - ts)
        tpos = te
        qpos = qe

    if qpos < len(a.query):
        parts += [{'pos': ttrue, 'type': 'I', 'length': len(a.query) - qpos}]
    
    return start+1, end+1, parts

def get_pileup_alignment_data(ref: str, seqs: list[str]) -> dict:
    """Get a local tile for use in a higlass-pileup plot.

    :param ref: The reference to align to
    :param seqs: The sequences to align to the reference
    """
    local_data = {
      "type": 'local-tiles',
      "tilesetInfo": {
        'min_pos': [0],
        'max_pos': [len(ref)],
        'max_width': len(ref),
        'tile_size': len(ref),
        'chromsizes': [['a', len(ref)]],
        'max_zoom': 0,
        'max_tile_width': 100000,
        'format': 'subs'
      },
      "tiles": {
        '0.0': [],
      }
    }
    
    for i,seq in enumerate(seqs):
        a = align_sequences(ref, seq)
        start, end, subs = get_subs(a)

        local_data['tiles']['0.0'].append({
          "id": f"r{i}",
          "from": start,
          "to": end,
          "substitutions": subs,
          "color": 0
        })

    return local_data