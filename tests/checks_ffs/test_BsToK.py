from slophep.FormFactors import BsToKFF
import check_utils as chk


def test_BCL_hammer():
    fflist = ["f0", "f+"]
    qsq, hammer_FFs = chk.read_txt_data("data/checks_hammer_BsToKFFBCL.txt", fflist)
    FF = BsToKFF.BCL()
    slop_FFs = chk.get_spectrum_slop(FF, qsq, "get_ff")

    for iff in fflist:
        assert chk.within_tolerance(hammer_FFs[iff], slop_FFs[iff])
