from clodius.tiles.pileup import align_sequences, get_subs, get_pileup_alignment_data

def test_get_subs():
    a = align_sequences("TTTTT", "AAAATTATTAAAA")
    s = get_subs(a)

    assert s[2][0]['type'] == 'I'
    assert s[2][0]['pos'] == 0
    assert s[2][0]['length'] == 4

    assert s[2][-1]['type'] == 'I'
    assert s[2][-1]['pos'] == 5
    assert s[2][-1]['length'] == 4

    a = align_sequences("TTTTT", "TTATT")
    s = get_subs(a)

    # assert 1-based start positions and closed intervals
    assert s[0] == 1
    assert s[1] == 6
    assert s[2][0]['pos'] == 2 # subs are 0-based
    assert s[2][0]['base'] == 'T'
    assert s[2][0]['variant'] == 'A'

    a = align_sequences("TTTTT", "TTATTT")
    s = get_subs(a)

    assert s[0] == 1
    assert s[1] == 6
    assert s[2][0]['pos'] == 2
    assert s[2][0]['type'] == 'I'
    assert s[2][0]['length'] == 1

def test_get_pileup_alignment_data():
    d = get_pileup_alignment_data(
        "ATTGA", ["TATTTTGGACCGCGCGTTCATTTACACGTC"]
    )
    
    assert 'type' in d
    assert 'tiles' in d
    assert d['tiles']['0.0'][0]['substitutions'][0]['type'] == 'I'