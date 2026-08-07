from slophep.FormFactors import BuToD1FF
import check_utils as chk


def test_BLR_hammer():
    fflist = ["fS", "fV1", "fV2", "fV3", "fA", "fT1", "fT2", "fT3"]
    qsq, hammer_FFs = chk.read_txt_data("data/checks_hammer_BuToD1FFBLR.txt", fflist)
    FF = BuToD1FF.BLR()
    slop_FFs = chk.get_spectrum_slop(FF, qsq, "get_ff")

    for iff in fflist:
        assert chk.within_tolerance(hammer_FFs[iff], slop_FFs[iff])
