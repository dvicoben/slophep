from math import sqrt

from slophep.FormFactors.FormFactorsBaryonic.FFOneHalfpToOneHalfpBase import FormFactorOneHalfpToOneHalfp
from slophep.Tools.errfluct_tools import FluctType, fluctsettings

class FFBaryonic_DKMR(FormFactorOneHalfpToOneHalfp):
    _name = "FFBaryonic@DKMR"

    def define_userparams(self):
        ffpar = {
            "a_Vt_0"  : 0.74392566003686,
            "a_Vt_1"  : -4.6476511382467,
            "a_Vt_2"  : 0.0,
            "a_At_0"  : 0.73960499809758,
            "a_At_1"  : -4.3664554520299,
            "a_At_2"  : 0.0,
            "a_Vl_0"  : 0.81458549275845,
            "a_Vl_1"  : -4.8987680713973,
            "a_Vl_2"  : 0.0,
            "a_Vp_0"  : 1.077996185173,
            "a_Vp_1"  : -6.4170836206134,
            "a_Vp_2"  : 0.0,
            "a_Al_0"  : 0.6846557018386,
            "a_Al_1"  : -4.4311552783827,
            "a_Al_2"  : 0.0,
            "a_Ap_1"  : -4.4633624514401,
            "a_Ap_2"  : 0.0,
            "a_Tl_0"  : 0.97518189850197,
            "a_Tl_1"  : -5.4999842229709,
            "a_Tl_2"  : 0.0,
            "a_Tp_0"  : 0.70539533108116,
            "a_Tp_1"  : -4.3577976693726,
            "a_Tp_2"  : 0.0,
            "a_T5l_0" : 0.67275080357035,
            "a_T5l_1" : -4.4321757873186,
            "a_T5l_2" : 0.0,
            "a_T5p_1" : -4.4927847697461,
            "a_T5p_2" : 0.0,
            # internalparams
            "m1p" : 6.768,
            "m1m" : 6.332,
            "m0p" : 6.725,
            "m0m" : 6.276
        }
        return ffpar

    def z(self, q2: float, mres: float):
        # eqns 3.3-3.5 in https://arxiv.org/pdf/1702.02243
        mB = self.get_param(f"m_{self.B}")
        mM = self.get_param(f"m_{self.C}")
        t0 = (mB-mM)**2
        tp = (mres)**2
        sq2 = sqrt(tp-q2)
        st0 = sqrt(tp-t0)
        return (sq2-st0)/(sq2+st0)

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff(self, q2: float) -> dict[str, float]:
        """DKMR form factors, https://arxiv.org/abs/1702.02243
        Implementation reproduces EOS https://github.com/eos/eos/blob/v1.0.13/eos/form-factors/parametric-dkmr2017-impl.hh

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict[str, float]
            The FFs
        """
        m0p = self.get_userparam("m0p")
        z0p = self.z(q2, m0p)
        m0m = self.get_userparam("m0m")
        z0m = self.z(q2, m0m)
        m1p = self.get_userparam("m1p")
        z1p = self.z(q2, m1p)
        m1m = self.get_userparam("m1m")
        z1m = self.z(q2, m1m)

        fV_t = 1.0/(1. - (q2/(m0p**2)))*(
            self.get_userparam("a_Vt_0")
            + self.get_userparam("a_Vt_1")*z0p
            + self.get_userparam("a_Vt_2")*(z0p**2)
        )
        fV_l = 1.0/(1. - (q2/(m1m**2)))*(
            self.get_userparam("a_Vl_0")
            + self.get_userparam("a_Vl_1")*z1m
            + self.get_userparam("a_Vl_2")*(z1m**2)
        )
        fV_p = 1.0/(1. - (q2/(m1m**2)))*(
            self.get_userparam("a_Vp_0")
            + self.get_userparam("a_Vp_1")*z1m
            + self.get_userparam("a_Vp_2")*(z1m**2)
        )
        fA_t = 1.0/(1. - (q2/(m0m**2)))*(
            self.get_userparam("a_At_0")
            + self.get_userparam("a_At_1")*z0m
            + self.get_userparam("a_At_2")*(z0m**2)
        )
        fA_l = 1.0/(1. - (q2/(m1p**2)))*(
            self.get_userparam("a_Al_0")
            + self.get_userparam("a_Al_1")*z1p
            + self.get_userparam("a_Al_2")*(z1p**2)
        )
        fA_p = 1.0/(1. - (q2/(m1p**2)))*(
            self.get_userparam("a_Al_0")
            + self.get_userparam("a_Ap_1")*z1p
            + self.get_userparam("a_Ap_2")*(z1p**2)
        )
        fT_l = 1.0/(1. - (q2/(m1m**2)))*(
            self.get_userparam("a_Tl_0")
            + self.get_userparam("a_Tl_1")*z1m
            + self.get_userparam("a_Tl_2")*(z1m**2)
        )
        fT_p = 1.0/(1. - (q2/(m1m**2)))*(
            self.get_userparam("a_Tp_0")
            + self.get_userparam("a_Tp_1")*z1m
            + self.get_userparam("a_Tp_2")*(z1m**2)
        )
        fT5_l = 1.0/(1. - (q2/(m1p**2)))*(
            self.get_userparam("a_T5l_0")
            + self.get_userparam("a_T5l_1")*z1p
            + self.get_userparam("a_T5l_2")*(z1p**2)
        )
        fT5_p = 1.0/(1. - (q2/(m1p**2)))*(
            self.get_userparam("a_T5l_0")
            + self.get_userparam("a_T5p_1")*z1p
            + self.get_userparam("a_T5p_2")*(z1p**2)
        )

        ff = {
            "Vt"  : fV_t,
            "V0"  : fV_l,
            "Vp"  : fV_p,
            "At"  : fA_t,
            "A0"  : fA_l,
            "Ap"  : fA_p,
            "T0"  : fT_l,
            "T50" : fT5_l,
            "Tp"  : fT_p,
            "T5p" : fT5_p
        }
        return ff