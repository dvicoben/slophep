import numpy as np
from slophep.FormFactors.FormFactorsBToP.FFBToPBase import FormFactorBToP
from slophep.Tools.SamplingTools import fluctsettings, FluctType
from slophep.FormFactors import hqet
import slophep.Math.derivatives as md
from flavio.physics.bdecays.formfactors.b_p.cln import h_to_f

import logging
logger = logging.getLogger(__name__)


def calc_Li1(w: float, IWs: dict[str, float]) -> list[float]:
    li1 = [
        0.0,
        4.*(3.* IWs["CHI3"] - (w-1.)*IWs["CHI2"]),
        -4.*IWs["CHI3"],
        4.*IWs["CHI2"],
        2.* IWs["ETA"] - 1.,
        -1.0,
        -2.* (IWs["ETA"] + 1.)/(w + 1.)
    ]
    return li1

def calc_Li2(w: float, IWs: dict[str, float], la2OverlaB2: float) -> list[float]:
    Li2 = [
        0.0,
        2.*IWs["BETA1"] + 4.*(3.*IWs["BETA3"] - (w-1.)*IWs["BETA2"]),
        2.*IWs["BETA1"] - 4.*IWs["BETA3"],
        4.*IWs["BETA2"],
        3*la2OverlaB2 + 2.*(w+1.)*IWs["PHI1"],
        la2OverlaB2 + 2.*(w+1.)*IWs["PHI1"],
        4.*IWs["PHI1"],
    ]
    return Li2

def calc_Mi(w: float, IWs: dict[str, float], la1OverlaB2: float, la2OverlaB2: float) -> list[float]:
    Mi = np.zeros(25) 
    Mi[8] = (
        la1OverlaB2+ 6.* la2OverlaB2/(w + 1.)
        -2.*(w-1.)*IWs["PHI1"] 
        - 2.*(2.*IWs["ETA"] - 1.)*(w - 1.)/(w + 1.)
    )
    Mi[9] = 3.*la2OverlaB2/(w+1.) + 2.*IWs["PHI1"] - (2.*IWs["ETA"] - 1.) * (w - 1.)/(w + 1.)
    Mi[10] = (
        la1OverlaB2/3. 
        - la2OverlaB2*(w+4.)/(2.*(w+1.))
        + 2.*(w + 2.)*IWs["PHI1Q"] - (2.*IWs["ETA"] - 1.)/(w + 1.)
    )
    return Mi


class FFBToP_BLPRXP(FormFactorBToP):
    _name = "FFBToP@BLPRXP"
    def define_userparams(self):
        ffpar = {
            "RhoStSq" : 1.104,
            "cSt"     : 2.392,
            "Chi21"   : -0.116,
            "Chi2p"   : 0.0,
            "Chi3p"   : 0.0,
            "Eta1"    : 0.336,
            "Etap"    : 0.0,
            "Phi1p"   : 0.252,
            "Beta21"  : 0.0,
            "Beta3p"  : 0.0,
            # internalparams
            "a"       : 1.509/np.sqrt(2.0),
            "as"      : 0.27,
            "as"      : 0.27,
            "mb"      : 4.707,
            "mc"      : 4.707 - 3.407,
            "mBBar"   : 5.313,
            "mDBar"   : 1.973,
            "rho1"    : -0.363,
            "la2"     : 0.12,
            "With1OverMb2"  : True,
            "With1OverMbMc" : True,
            "WithAsOverM"   : True
        }
        return ffpar

    def get_internal(self) -> dict[str, float]:
        _as    = self.get_userparam("as")
        ash   = _as/np.pi
        mb    = self.get_userparam("mb")
        mc    = self.get_userparam("mc")
        mBBar = self.get_userparam("mBBar")
        mDBar = self.get_userparam("mDBar")
        rho1  = self.get_userparam("rho1")
        la2   = self.get_userparam("la2")
        # Some calculations       
        zBC = mc/mb
        laB = (mb*mBBar - mc*mDBar)/(mb-mc) - (mb+mc) + rho1/(4.*mb*mc)
        la1 = (2.*mb*mc)/(mb - mc)*(mBBar - mDBar - (mb-mc)) + rho1*(mb+mc)/(2.*mb*mc)
        corr1S = 2.* np.power(_as/3. * 0.796, 2.)
        dmbc = mb - mc
        mb_1 = mb * corr1S
        mc_1 = mb * corr1S
        LambdaBar_0 = (mb*mBBar - mc*mDBar)/dmbc - (2.*mb-dmbc) + rho1/(4*mb*mc)
        LambdaBar_1 = (mb_1 * mBBar - mc_1 * mDBar)/dmbc - 2.*mb_1
        upsilonb = LambdaBar_1/LambdaBar_0 - mb_1/mb
        upsilonc = LambdaBar_1/LambdaBar_0 - mc_1/mc

        params = {
            "ash"      : ash,
            "zBC"      : zBC,
            "cmagb"    : 3./2.*(0.5*np.log(zBC) + 13./9.),
            "cmagc"    : -3./2.*(0.5*np.log(zBC) + 13./9.),
            "laB"      : laB,
            "la1"      : la1,
            "la2/laB2" : la2/(laB*laB),
            "la1/laB2" : la1/(laB*laB),
            "corr1S"   : corr1S,
            "dmbc"     : dmbc,
            "mb1"      : mb_1,
            "mc1"      : mc_1,
            "LamBar0"  : LambdaBar_0,
            "LamBar1"  : LambdaBar_1,
            "eb"       : laB/(2.*mb),
            "ec"       : laB/(2.*mc),
            "upsilonb" : upsilonb,
            "upsilonc" : upsilonc
        }
        return params

    def z_of_w(self, w: float, a: float) -> float:
        sqrt2 = np.sqrt(2)
        return (np.sqrt(w+1) - sqrt2*a)/(np.sqrt(w+1) + sqrt2*a)

    @fluctsettings(FluctType.NUMERIC)
    def xi(self, w: float) -> float:
        a = self.get_userparam("a")
        a2 = a**2
        a4 = a**4
        zCon = self.z_of_w(w, a)
        zCon1 = (1.-a)/(1.+a)
        RhoSq = self.get_userparam("RhoStSq")
        cSt = self.get_userparam("cSt")
        Xi = (
            (1. - 8*a2*RhoSq*zCon + 16.*(2*cSt * a4 - RhoSq * a2)*zCon*zCon)
            /(1. - 8*a2*RhoSq*zCon1 + 16.*(2*cSt * a4 - RhoSq * a2)*zCon1*zCon1)
        )
        return Xi

    def IW(self, w: float, auxparams: dict[str, float]) -> dict:
        IWs = {
            "CHI2"  : self.get_userparam("Chi21") + (w-1.0)*self.get_userparam("Chi2p"),
            "CHI3"  : self.get_userparam("Chi3p")*(w-1.0),
            "ETA"   : self.get_userparam("Eta1") + (w-1.0)*self.get_userparam("Etap"),
            "BETA1" : auxparams["la1/laB2"]/4.0,
            "BETA2" : self.get_userparam("Beta21"),
            "BETA3" : auxparams["la2/laB2"]/8.0 + self.get_userparam("Beta3p")*(w-1.0),
            "PHI1"  : (auxparams["la1/laB2"]/3. - auxparams["la2/laB2"]/2.)/2. + (self.get_userparam("Phi1p"))*(w-1.),
            "PHI1Q" : self.get_userparam("Phi1p")
        }
        return IWs

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_hhat(self, q2: float) -> dict[str, float]:
        auxparams = self.get_internal()
        w = max(self.w(q2), 1)
        zBC = auxparams["zBC"]

        # Hs = 1.
        Hp = 1.
        Hm = 0.
        Ht = 1.

        # QCD correction functions
        # Cps = hqet.CP(w, zBC)
        Cv1 = hqet.CV1(w, zBC)
        Cv2 = hqet.CV2(w, zBC)
        Cv3 = hqet.CV3(w, zBC)
        Ct1 = hqet.CT1(w, zBC)
        Ct2 = hqet.CT2(w, zBC)
        Ct3 = hqet.CT3(w, zBC)

        # alpha_s
        ash = auxparams["ash"]
        Hp += ash*(Cv1+0.5*(w+1)*(Cv2+Cv3))
        Hm += ash*0.5*(w+1)*(Cv2-Cv3)
        Ht += ash*(Ct1-Ct2+Ct3)

        # Epsilon b and epsilon c
        eb = auxparams["eb"]
        ec = auxparams["ec"]

        IWs = self.IW(w)
        Li1 = calc_Li1(w, IWs)
        Hp += (ec + eb)*Li1[1]
        Hm += (ec - eb)*Li1[4]
        Ht += (ec + eb)*(Li1[1]-Li1[4])

        # Upsilon Expansion
        corrb = eb*auxparams["upsilonb"]
        corrc = ec*auxparams["upsilonc"]
        Hp += 0.
        Hm += -(corrc-corrb)
        Ht += (corrc+corrb)

        # Epsilon c squared
        ec2 = ec*ec
        Li2 = calc_Li2(w, IWs, auxparams["la2/laB2"])
        # Hps += ec2*(Li2[2]+Li2[3]*(w-1)+Li2[5]-Li2[6]*(w+1))
        Hp += ec2*Li2[1]
        Hm += ec2*Li2[4]
        Ht += ec2*(Li2[1]-Li2[4])


        # Epsilon b squared
        if self.get_userparam("With1OverMb2"):
            eb2 = eb*eb
            Hp += +eb2*Li2[1]
            Hm += -eb2*Li2[4]
            Ht += +eb2*(Li2[1]-Li2[4])
        # Epsilon b * Epsilon c
        if self.get_userparam("With1OverMbMc"):
            eceb = ec*eb
            Mi = calc_Mi(w, IWs)
            Hp += eceb * (Mi[1] - Mi[8])
            Hm += 0.
            Ht += eceb * (Mi[1] + Mi[8] - 2.*Mi[15])
        # alpha_s/mb, alpha_s/mc
        if self.get_userparam("WithAsOverM"):
            asec = ash*ec
            aseb = ash*eb
            cmagc = auxparams["cmagc"]
            cmagb = auxparams["cmagb"]
            bound_lo = 1.0 + md._DER_EPS
            dCv1 = md.derivative(lambda x: hqet.CV1(x, zBC), w, bound_lo = bound_lo)
            dCv2 = md.derivative(lambda x: hqet.CV2(x, zBC), w, bound_lo = bound_lo)
            dCv3 = md.derivative(lambda x: hqet.CV3(x, zBC), w, bound_lo = bound_lo)
            dCt1 = md.derivative(lambda x: hqet.CT1(x, zBC), w, bound_lo = bound_lo)
            dCt2 = md.derivative(lambda x: hqet.CT2(x, zBC), w, bound_lo = bound_lo)
            dCt3 = md.derivative(lambda x: hqet.CT3(x, zBC), w, bound_lo = bound_lo)

            Hp += asec * (cmagc*Li1[1] + Li1[1]*Cv1 - Li1[5]*(-1 + w)*Cv3 + ((1 + w)*(Li1[1] - (Li1[4]*(-1 + w))/(1 + w))*(Cv2 + Cv3))/2. + 2*(-1 + w)*(dCv1 + ((1 + w)*(dCv2 + dCv3))/2.))
            Hm += asec * (Li1[4]*Cv1 + Li1[5]*(1 + w)*Cv3 + ((1 + w)*((Li1[1] - (Li1[4]*(-1 + w))/(1 + w))*(Cv2 - Cv3) + 2*(-1 + w)*(dCv2 - dCv3)))/2.)
            Ht += asec * (cmagc*Li1[1] - Li1[4]*Ct1 - 2*Li1[5]*Ct3 + Li1[1]*(Ct1 - Ct2 + Ct3) + Li1[4]*(Ct2 + Ct3) + 2*(-1 + w)*(dCt1 - dCt2 + dCt3))
            
            Hp += aseb * (cmagb*Li1[1] + Li1[1]*Cv1 - Li1[5]*(-1 + w)*Cv2 + ((1 + w)*(Li1[1] - (Li1[4]*(-1 + w))/(1 + w))*(Cv2 + Cv3))/2. + 2*(-1 + w)*(dCv1 + ((1 + w)*(dCv2 + dCv3))/2.))
            Hm += aseb * (-(Li1[4]*Cv1) - Li1[5]*(1 + w)*Cv2 + ((1 + w)*((Li1[1] - (Li1[4]*(-1 + w))/(1 + w))*(Cv2 - Cv3) + 2*(-1 + w)*(dCv2 - dCv3)))/2.)
            Ht += aseb * (cmagb*Li1[1] - Li1[4]*Ct1 + 2*Li1[5]*Ct2 + Li1[1]*(Ct1 - Ct2 + Ct3) - Li1[4]*(Ct2 + Ct3) + 2*(-1 + w)*(dCt1 - dCt2 + dCt3))

        h = {
            "h+" : Hp,
            "h-" : Hm,
            "hT" : Ht
        }
        return h

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff(self, q2: float) -> dict:
        """FF in BLPRXP parameterisation from https://arxiv.org/abs/2206.11281 as in HAMMER v1.4.1

        Parameters
        ----------
        q2 : float
            q2 value to calculate FF at

        Returns
        -------
        dict
            FF dictionary
        """
        Mb = self.get_param(f"m_{self.B}")
        Mc = self.get_param(f"m_{self.P}")
        w = max(self.w(q2), 1)
        # IG
        Xi = self.xi(w)
        hhat = self.get_hhat(q2)
        h = {iff : Xi*hhat[iff] for iff in hhat}
        return h_to_f(Mb, Mc, h, q2)


    