import numpy as np
from slophep.FormFactors.FormFactorsBToDstst.FFBToDststBase import (
    FormFactorBToD0st, FormFactorBToD1, FormFactorBToD1st, FormFactorBToD2st
)
from slophep.FormFactors import hqet
from slophep.Tools.errfluct_tools import fluctsettings, FluctType

import logging
logger = logging.getLogger(__name__)


class FFBToD1_BLR(FormFactorBToD1):
    _name = "FFBToD1@BLR"

    def define_userparams(self):
        ffpar = {
            "t1"   : 0.7 ,
            "tp"   : -1.6,
            "tau1" : -0.5,
            "tau2" : 2.9 ,
            "eta1" : 0.  ,
            "eta2" : 0.  ,
            "eta3" : 0.  ,
            #internalparams
            "ash"  : 0.26/np.pi   ,
            "mb"   : 4.710        ,
            "mc"   : 4.710 - 3.400,
            "laB"  : 0.4          ,
            "laP"  : 0.8          ,
        }
        return ffpar

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff_mmeson(self, q2: float, mC: float, mB: float = None) -> dict:
        """
        BLR form factors as in https://arxiv.org/pdf/1711.03110,
        https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoD1BLR.cc
        """
        Mb = mB if mB is not None else self.get_param(f"m_{self.B}")
        Mc = mC

        Mb2 = Mb*Mb
        Mc2 = Mc*Mc
        w = max((Mb2 + Mc2 - q2)/(2.*Mb*Mc), 1.0)

        eB   = 0.5 / self.get_userparam("mb")
        eC   = 0.5 / self.get_userparam("mc")
        zBC = self.get_userparam("mc")/self.get_userparam("mb")
        ash  = self.get_userparam("ash")
        laB  = self.get_userparam("laB")
        laP  = self.get_userparam("laP")
        t1   = self.get_userparam("t1")
        tp   = self.get_userparam("tp")
        tau1 = self.get_userparam("tau1")
        tau2 = self.get_userparam("tau2")
        eta1 = self.get_userparam("eta1")
        eta2 = self.get_userparam("eta2")
        eta3 = self.get_userparam("eta3")

        LambdaD32 = -laB + laP*w
        Fb = laB + laP - tau2 - tau1*(1 + 2*w)
        LOIWtau = t1 + (w-1)*t1*tp

        # QCD correction functions
        Cs  = hqet.CS(w, zBC)
        Cv1 = hqet.CV1(w, zBC)
        Cv2 = hqet.CV2(w, zBC)
        Cv3 = hqet.CV3(w, zBC)
        Ca1 = hqet.CA1(w, zBC)
        Ct1 = hqet.CT1(w, zBC)
        Ct2 = hqet.CT2(w, zBC)
        Ct3 = hqet.CT3(w, zBC)

        Fsc = (-(eC*(4*LambdaD32 + 2*(1 + w)*(6*eta1 + 2*(w-1)*eta2 - eta3) - 2*(w-1)*((1 + 2*w)*tau1 + tau2))) - 2*(1 + w)*(1 + ash*Cs) - 2*(w-1)*eB*Fb)/np.sqrt(6.)
        Fv1 = (-(eC*(4*(1 + w)*LambdaD32 - (-1 + w*w)*(2*eta1 + 3*eta3 + 3*tau1 - 3*tau2))) + (1 - w*w)*(1 + ash*Cv1) - (-1 + w*w)*eB*Fb)/np.sqrt(6.)
        Fv2 = (-3 - eC*(10*eta1 + 4*(w-1)*eta2 - 5*eta3 + (-1 + 4*w)*tau1 + 5*tau2) - ash*(3*Cv1 + 2*(1 + w)*Cv2) - 3*eB*Fb)/np.sqrt(6.)
        Fv3 = (-2 + w + eC*(4*LambdaD32 - 2*(6 + w)*eta1 - 4*(w-1)*eta2 - (-2 + 3*w)*eta3 + (2 + w)*tau1 + (2 + 3*w)*tau2) - ash*((2 - w)*Cv1 + 2*(1 + w)*Cv3) + (2 + w)*eB*Fb)/np.sqrt(6.)
        Fa = (-(eC*(4*LambdaD32 - (1 + w)*(2*eta1 + 3*eta3) - 3*(w-1)*(tau1 - tau2))) + (-1 - w)*(1 + ash*Ca1) - (w-1)*eB*Fb)/np.sqrt(6.)
        Ft1 = (-(eC*(4*LambdaD32 + (1 + w)*(2*eta1 + 3*eta3) - 3*(w-1)*(tau1 - tau2))) + (1 + w)*(1 + ash*(Ct1 + (w-1)*Ct2)) + (w-1)*eB*Fb)/np.sqrt(6.)
        Ft2 = (-(eC*(4*LambdaD32 - (1 + w)*(2*eta1 + 3*eta3) - 3*(w-1)*(tau1 - tau2))) + (-1 - w)*(1 + ash*(Ct1 - (w-1)*Ct3)) + (w-1)*eB*Fb)/np.sqrt(6.)
        Ft3 = (3 - eC*(-10*eta1 - 4*(w-1)*eta2 + 5*eta3 + (-1 + 4*w)*tau1 + 5*tau2) + ash*(3*Ct1 - (2 - w)*Ct2 + 3*Ct3) + 3*eB*Fb)/np.sqrt(6.)

        ffs = {
            "fS"  : LOIWtau*Fsc,
            "fV1" : LOIWtau*Fv1,
            "fV2" : LOIWtau*Fv2,
            "fV3" : LOIWtau*Fv3,
            "fA"  : LOIWtau*Fa ,
            "fT1" : LOIWtau*Ft1,
            "fT2" : LOIWtau*Ft2,
            "fT3" : LOIWtau*Ft3,
        }
        return ffs



class FFBToD2st_BLR(FormFactorBToD2st):
    _name = "FFBToD2st@BLR"

    def define_userparams(self):
        ffpar = {
            "t1"   : 0.7 ,
            "tp"   : -1.6,
            "tau1" : -0.5,
            "tau2" : 2.9 ,
            "eta1" : 0.  ,
            "eta2" : 0.  ,
            "eta3" : 0.  ,
            #internalparams
            "ash"  : 0.26/np.pi   ,
            "mb"   : 4.710        ,
            "mc"   : 4.710 - 3.400,
            "laB"  : 0.4          ,
            "laP"  : 0.8          ,
        }
        return ffpar

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff_mmeson(self, q2: float, mC: float, mB: float = None) -> dict:
        """
        BLR form factors as in https://arxiv.org/pdf/1711.03110,
        https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoD1BLR.cc
        """
        Mb = mB if mB is not None else self.get_param(f"m_{self.B}")
        Mc = mC

        Mb2 = Mb*Mb
        Mc2 = Mc*Mc
        w = max((Mb2 + Mc2 - q2)/(2.*Mb*Mc), 1.0)

        eB   = 0.5 / self.get_userparam("mb")
        eC   = 0.5 / self.get_userparam("mc")
        zBC = self.get_userparam("mc")/self.get_userparam("mb")
        ash  = self.get_userparam("ash")
        laB  = self.get_userparam("laB")
        laP  = self.get_userparam("laP")
        t1   = self.get_userparam("t1")
        tp   = self.get_userparam("tp")
        tau1 = self.get_userparam("tau1")
        tau2 = self.get_userparam("tau2")
        eta1 = self.get_userparam("eta1")
        eta2 = self.get_userparam("eta2")
        eta3 = self.get_userparam("eta3")

        Fb = laB + laP - tau2 - tau1*(1 + 2*w)
        LOIWtau = t1 + (w-1)*t1*tp

        # QCD correction functions
        Cps = hqet.CP(w, zBC)
        Cv1 = hqet.CV1(w, zBC)
        Ca1 = hqet.CA1(w, zBC)
        Ca2 = hqet.CA2(w, zBC)
        Ca3 = hqet.CA3(w, zBC)
        Ct1 = hqet.CT1(w, zBC)
        Ct2 = hqet.CT2(w, zBC)
        Ct3 = hqet.CT3(w, zBC)

        Kp = 1 + eC*(-2*eta1 - 2*(w-1)*eta2 + eta3 + (1 + 2*w)*tau1 + tau2) + ash*Cps + eB*Fb
        Ka1 = -(eC*(-((1 + w)*(2*eta1 - eta3)) + (w-1)*(tau1 - tau2))) + (-1 - w)*(1 + ash*Ca1) - (w-1)*eB*Fb
        Ka2 = -2*eC*(eta2 + tau1) + ash*Ca2
        Ka3 = 1 - eC*(2*eta1 - 2*eta2 - eta3 + tau1 + tau2) + ash*(Ca1 + Ca3) + eB*Fb
        Kv = -1 - eC*(-2*eta1 + eta3 + tau1 - tau2) - ash*Cv1 - eB*Fb
        Kt1 = 1 - eC*(2*eta1 - eta3) + ash*(Ct1 + ((w-1)*(Ct2 - Ct3))/2.)
        Kt2 = -(eC*(tau1 - tau2)) + ((1 + w)*ash*(Ct2 + Ct3))/2. + eB*Fb
        Kt3 = 2*eC*(-eta2 + tau1) - ash*Ct2

        ffs = {
            "kP"  : LOIWtau*Kp ,
            "kA1" : LOIWtau*Ka1,
            "kA2" : LOIWtau*Ka2,
            "kA3" : LOIWtau*Ka3,
            "kV"  : LOIWtau*Kv ,
            "kT1" : LOIWtau*Kt1,
            "kT2" : LOIWtau*Kt2,
            "kT3" : LOIWtau*Kt3,
        }
        return ffs



class FFBToD0st_BLR(FormFactorBToD0st):
    _name = "FFBToD0st@BLR"

    def define_userparams(self):
        ffpar = {
            "zt1"   : 0.7,
            "ztp"   : 0.2,
            "zeta1" : 0.6,
            "chi1"  : 0.0,
            "chi2"  : 0.0,
            #internalparams
            "ash" : 0.26/np.pi   ,
            "mb"  : 4.710        ,
            "mc"  : 4.710 - 3.400,
            "laB" : 0.4          ,
            "laS" : 0.76         ,
        }
        return ffpar

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff_mmeson(self, q2: float, mC: float, mB: float = None) -> dict:
        """
        BLR form factors as in https://arxiv.org/pdf/1711.03110,
        https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoD1BLR.cc
        """
        Mb = mB if mB is not None else self.get_param(f"m_{self.B}")
        Mc = mC

        Mb2 = Mb*Mb
        Mc2 = Mc*Mc
        w = max((Mb2 + Mc2 - q2)/(2.*Mb*Mc), 1.0)

        mb    = self.get_userparam("mb")
        mc    = self.get_userparam("mc")
        eB    = 0.5 / mb
        eC    = 0.5 / mc
        zBC   = mc/mb
        ash   = self.get_userparam("ash")
        laB   = self.get_userparam("laB")
        laS   = self.get_userparam("laS")
        zt1   = self.get_userparam("zt1")
        ztp   = self.get_userparam("ztp")
        zeta1 = self.get_userparam("zeta1")
        chi1  = self.get_userparam("chi1")
        chi2  = self.get_userparam("chi2")

        LambdaD12 = -laB + laS*w
        Gb = (-(laB*(2 + w)) + laS*(1 + 2*w))/(1 + w) - 2*(w-1)*zeta1
        LOIWzeta = zt1 + (w-1)*zt1*ztp

        # QCD correction functions
        Cps = hqet.CP(w, zBC)
        Ca1 = hqet.CA1(w, zBC)
        Ca2 = hqet.CA2(w, zBC)
        Ca3 = hqet.CA3(w, zBC)
        Ct1 = hqet.CT1(w, zBC)
        
        gps = eC*(3*LambdaD12 - 2*(-1 + w*w)*zeta1 + (w-1)*(6*chi1 - 2*(1 + w)*chi2)) + (w-1)*(1 + ash*Cps) - (1 + w)*eB*Gb
        gp = -(eC*((3*LambdaD12)/(1 + w) - 2*(w-1)*zeta1)) + ((w-1)*ash*(Ca2 + Ca3))/2. - eB*Gb
        gm = 1 + eC*(6*chi1 - 2*(1 + w)*chi2) + ash*(Ca1 + ((w-1)*(Ca2 - Ca3))/2.)
        gt = 1 + eC*((3*LambdaD12)/(1 + w) - 2*(w-1)*zeta1 + 6*chi1 - 2*(1 + w)*chi2) + ash*Ct1 - eB*Gb

        ffs = {
            "gP" : LOIWzeta*gps,
            "g+" : LOIWzeta*gp,
            "g-" : LOIWzeta*gm,
            "gT" : LOIWzeta*gt,
        }
        return ffs



class FFBToD1st_BLR(FormFactorBToD1st):
    _name = "FFBToD1st@BLR"

    def define_userparams(self):
        ffpar = {
            "zt1"   : 0.7,
            "ztp"   : 0.2,
            "zeta1" : 0.6,
            "chi1"  : 0.0,
            "chi2"  : 0.0,
            #internalparams
            "ash" : 0.26/np.pi   ,
            "mb"  : 4.710        ,
            "mc"  : 4.710 - 3.400,
            "laB" : 0.4          ,
            "laS" : 0.76         ,
        }
        return ffpar

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff_mmeson(self, q2: float, mC: float, mB: float = None) -> dict:
        """
        BLR form factors as in https://arxiv.org/pdf/1711.03110,
        https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoD1BLR.cc
        """
        Mb = mB if mB is not None else self.get_param(f"m_{self.B}")
        Mc = mC

        Mb2 = Mb*Mb
        Mc2 = Mc*Mc
        w = max((Mb2 + Mc2 - q2)/(2.*Mb*Mc), 1.0)

        mb    = self.get_userparam("mb")
        mc    = self.get_userparam("mc")
        eB    = 0.5 / mb
        eC    = 0.5 / mc
        zBC   = mc/mb
        ash   = self.get_userparam("ash")
        laB   = self.get_userparam("laB")
        laS   = self.get_userparam("laS")
        zt1   = self.get_userparam("zt1")
        ztp   = self.get_userparam("ztp")
        zeta1 = self.get_userparam("zeta1")
        chi1  = self.get_userparam("chi1")
        chi2  = self.get_userparam("chi2")

        LambdaD12 = -laB + laS*w
        Gb = (-(laB*(2 + w)) + laS*(1 + 2*w))/(1 + w) - 2*(w-1)*zeta1
        LOIWzeta = zt1 + (w-1)*zt1*ztp

        # QCD correction functions
        Cs  = hqet.CS(w, zBC)
        Cv1 = hqet.CV1(w, zBC)
        Cv2 = hqet.CV2(w, zBC)
        Cv3 = hqet.CV3(w, zBC)
        Ca1 = hqet.CA1(w, zBC)
        Ct1 = hqet.CT1(w, zBC)
        Ct2 = hqet.CT2(w, zBC)
        Ct3 = hqet.CT3(w, zBC)

        Gs  = 1 - eC*(LambdaD12/(1 + w) - 2*(w-1)*zeta1 + 2*chi1 - 2*(1 + w)*chi2) + ash*Cs - eB*Gb
        Gv1 = eC*(LambdaD12 - 2*(w-1)*chi1) + (w-1)*(1 + ash*Cv1) - (1 + w)*eB*Gb
        Gv2 = eC*(2*zeta1 - 2*chi2) - ash*Cv2
        Gv3 = -1 - eC*(LambdaD12/(1 + w) + 2*zeta1 - 2*chi1 + 2*chi2) - ash*(Cv1 + Cv3) + eB*Gb
        Ga  = 1 + eC*(LambdaD12/(1 + w) - 2*chi1) + ash*Ca1 - eB*Gb
        Gt1 = -1 + eC*(LambdaD12/(1 + w) + 2*chi1) - ash*(Ct1 + (w-1)*Ct2) + eB*Gb
        Gt2 = 1 + eC*(LambdaD12/(1 + w) - 2*chi1) + ash*(Ct1 - (w-1)*Ct3) + eB*Gb
        Gt3 = eC*(2*zeta1 + 2*chi2) - ash*Ct2

        ffs = {
            "gS"  : LOIWzeta*Gs ,
            "gV1" : LOIWzeta*Gv1,
            "gV2" : LOIWzeta*Gv2,
            "gV3" : LOIWzeta*Gv3,
            "gA"  : LOIWzeta*Ga ,
            "gT1" : LOIWzeta*Gt1,
            "gT2" : LOIWzeta*Gt2,
            "gT3" : LOIWzeta*Gt3,
        }
        return ffs