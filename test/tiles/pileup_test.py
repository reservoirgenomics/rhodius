from clodius.tiles.pileup import get_pileup_alignment_data
from clodius.alignment import align_sequences, alignment_to_subs, order_by_clustering


def test_alignment_to_subs():
    a = align_sequences("TTTTT", "AAAATTATTAAAA")
    print("")
    print(a)
    s = alignment_to_subs(a)

    print("s", s)

    assert s[2][0]["type"] == "I"
    assert s[2][0]["pos"] == 0
    assert s[2][0]["length"] == 4

    assert s[2][-1]["type"] == "I"
    assert s[2][-1]["pos"] == 5
    assert s[2][-1]["length"] == 4

    a = align_sequences("TTTTT", "TTATT")
    s = alignment_to_subs(a)

    # assert 1-based start positions and closed intervals
    assert s[0] == 1
    assert s[1] == 6
    assert s[2][0]["pos"] == 2  # subs are 0-based
    assert s[2][0]["base"] == "T"
    assert s[2][0]["variant"] == "A"

    a = align_sequences("TTTTT", "TTATTT")
    s = alignment_to_subs(a)

    assert s[0] == 1
    assert s[1] == 6
    assert s[2][0]["pos"] == 2
    assert s[2][0]["type"] == "I"
    assert s[2][0]["length"] == 1
