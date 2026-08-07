from slophep.FormFactors import BdToDstFF
import check_utils as chk


def test_CLN_hammer():
    fflist = ["f", "g", "F1", "F2"]
    qsq, hammer_FFs = chk.read_txt_data("data/checks_hammer_BdToDstFFCLN.txt", ["Ff", "Fg", "Fm", "Fp"])
    FF = BdToDstFF.CLN()
    slop_FFs = chk.get_spectrum_slop(FF, qsq, "get_ff_gfF1F2_basis")

    # Need to convert hammer_FFs to common basis
    hammer_conv_FFs = chk.convert_btov_hammer_bglbasis(hammer_FFs, FF.get_param("m_B0"), FF.get_param("m_D*+"), qsq)

    for iff in fflist:
        assert chk.within_tolerance(hammer_conv_FFs[iff], slop_FFs[iff])


def test_BLPR_hammer():
    fflist = ["hA1", "hA2", "hA3", "hV", "hT1", "hT2", "hT3"]
    qsq, hammer_FFs = chk.read_txt_data("data/checks_hammer_BdToDstFFBLPR.txt", ["Fs", "Ff", "Fg", "Fm", "Fp", "Fzt", "Fmt", "Fpt"])
    FF = BdToDstFF.BLPR_Hammer()
    slop_FFs = chk.get_spectrum_slop(FF, qsq, "get_ff_h_basis")

    # Need to convert hammer_FFs to common basis
    hammer_conv_FFs = chk.convert_btov_hammer_hbasis(hammer_FFs, FF.get_param("m_B0"), FF.get_param("m_D*+"), qsq)

    for iff in fflist:
        assert chk.within_tolerance(hammer_conv_FFs[iff], slop_FFs[iff])


def test_BGL_hammer():
    fflist = ["f", "g", "F1", "F2"]
    qsq, hammer_FFs = chk.read_txt_data("data/checks_hammer_BdToDstFFBGL.txt", ["Ff", "Fg", "Fm", "Fp"])
    FF = BdToDstFF.BGL_Hammer()
    slop_FFs = chk.get_spectrum_slop(FF, qsq, "get_ff_gfF1F2_basis")

    # Need to convert hammer_FFs to common basis
    hammer_conv_FFs = chk.convert_btov_hammer_bglbasis(hammer_FFs, FF.get_param("m_B0"), FF.get_param("m_D*+"), qsq)

    for iff in fflist:
        assert chk.within_tolerance(hammer_conv_FFs[iff], slop_FFs[iff])



def test_BSZ_eos():
    fflist = ["A0", "A1", "A12", "V", "T1", "T2", "T23"]
    qsq, eos_FFs = chk.read_txt_data("data/checks_eos_BdToDstFFBSZ.txt", fflist)
    FF = BdToDstFF.BSZ()
    slop_FFs = chk.get_spectrum_slop(FF, qsq, "get_ff")
    
    for iff in fflist:
        # Increase tolerance a bit (to ~0.05% from 0.01%) to accommodate differences induced by slightly different masses
        assert chk.within_tolerance(eos_FFs[iff], slop_FFs[iff], tol=5e-4)