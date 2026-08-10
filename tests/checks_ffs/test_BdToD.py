from slophep.FormFactors import BdToDFF
import check_utils as chk


def test_BLPR_hammer():
    fflist = ["f0", "f+", "fT"]
    qsq, hammer_FFs = chk.read_txt_data("data/checks_hammer_BdToDFFBLPR.txt", ["Fs", "f0", "f+", "fT"])
    FF = BdToDFF.BLPR_Hammer()
    slop_FFs = chk.get_spectrum_slop(FF, qsq, "get_ff")

    for iff in fflist:
        assert chk.within_tolerance(hammer_FFs[iff], slop_FFs[iff])