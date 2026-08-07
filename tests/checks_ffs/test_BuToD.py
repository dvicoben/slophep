from slophep.FormFactors import BuToDFF
import check_utils as chk


def test_CLN_hammer():
    fflist = ["f0", "f+"]
    qsq, hammer_FFs = chk.read_txt_data("data/checks_hammer_BuToDFFCLN.txt", ["Fs", "f0", "f+"])
    FF = BuToDFF.CLN()
    slop_FFs = chk.get_spectrum_slop(FF, qsq, "get_ff")

    for iff in fflist:
        assert chk.within_tolerance(hammer_FFs[iff], slop_FFs[iff])


def test_BGL_hammer():
    fflist = ["f0", "f+"]
    qsq, hammer_FFs = chk.read_txt_data("data/checks_hammer_BuToDFFBGL.txt", ["Fs", "f0", "f+"])
    FF = BuToDFF.BGL_Hammer()
    slop_FFs = chk.get_spectrum_slop(FF, qsq, "get_ff")

    for iff in fflist:
        assert chk.within_tolerance(hammer_FFs[iff], slop_FFs[iff])


def test_BLPR_hammer():
    fflist = ["f0", "f+", "fT"]
    qsq, hammer_FFs = chk.read_txt_data("data/checks_hammer_BuToDFFBLPR.txt", ["Fs", "f0", "f+", "fT"])
    FF = BuToDFF.BLPR_Hammer()
    slop_FFs = chk.get_spectrum_slop(FF, qsq, "get_ff")

    for iff in fflist:
        assert chk.within_tolerance(hammer_FFs[iff], slop_FFs[iff])


def test_BSZ_eos():
    fflist = ["f+", "f0", "fT"]
    qsq, eos_FFs = chk.read_txt_data("data/checks_eos_BuToDFFBSZ.txt", ["f0", "f+", "fT"])
    FF = BuToDFF.BSZ()
    slop_FFs = chk.get_spectrum_slop(FF, qsq, "get_ff")

    for iff in fflist:
        # Increase tolerance a bit (to ~0.05% from 0.01%) to accommodate differences induced by slightly different masses
        assert chk.within_tolerance(eos_FFs[iff], slop_FFs[iff], tol=5e-4)