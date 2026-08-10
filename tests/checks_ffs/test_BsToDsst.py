from slophep.FormFactors import BsToDsstFF
import check_utils as chk


def test_BGL_hammer():
    fflist = ["f", "g", "F1", "F2"]
    qsq, hammer_FFs = chk.read_txt_data("data/checks_hammer_BsToDsstFFBGL.txt", ["Ff", "Fg", "Fm", "Fp"])
    FF = BsToDsstFF.BGL_Hammer()
    slop_FFs = chk.get_spectrum_slop(FF, qsq, "get_ff_gfF1F2_basis")

    # Need to convert hammer_FFs to common basis
    hammer_conv_FFs = chk.convert_btov_hammer_bglbasis(hammer_FFs, FF.get_param("m_Bs"), FF.get_param("m_Ds*"), qsq)

    for iff in fflist:
        assert chk.within_tolerance(hammer_conv_FFs[iff], slop_FFs[iff])
