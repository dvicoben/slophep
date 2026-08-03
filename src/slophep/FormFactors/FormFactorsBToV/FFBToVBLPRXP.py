import numpy as np
from slophep.FormFactors.FormFactorsBToV.FFBToVBase import FormFactorBToV
from slophep.Tools.errfluct_tools import fluctsettings, FluctType
from slophep.FormFactors import hqet
import slophep.Math.derivatives as md
from flavio.physics.bdecays.formfactors.b_v.cln import h_to_A

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


class FFBToV_BLPRXP(FormFactorBToV):
    _name = "FFBToV@BLPRXP"
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
            "a"       : 1.509/np.sqrt(2.0),
            "as"      : 0.27,
            "mb"      : 4.707,
            "mc"      : 4.707 - 3.407,
            "mBBar"   : 5.313,
            "mDBar"   : 1.973,
            "rho1"    : -0.363,
            "la2"     : 0.12,
            "With1OverMb2"  : True,
            "With1OverMbMc" : True,
            "WithAsOverM"   : True,
            "WithHA1As2"    : True
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

    @fluctsettings(FluctType.DICTNUMERIC)
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
    def get_hhat(self, q2: float) -> dict:
        auxparams = self.get_internal()
        w = max(self.w(q2), 1)
        zBC = auxparams["zBC"]

        # Hps = 1.
        Hv  = 1.
        Ha1 = 1.
        Ha2 = 0.
        Ha3 = 1.
        Ht1 = 1.
        Ht2 = 0.
        Ht3 = 0.

        # QCD correction functions
        # Cps = hqet.CP(w, zBC)
        Cv1 = hqet.CV1(w, zBC)
        Ca1 = hqet.CA1(w, zBC)
        Ca2 = hqet.CA2(w, zBC)
        Ca3 = hqet.CA3(w, zBC)
        Ct1 = hqet.CT1(w, zBC)
        Ct2 = hqet.CT2(w, zBC)
        Ct3 = hqet.CT3(w, zBC)
        
        # alpha_s
        ash = auxparams["ash"]
        # Hps += ash*Cps
        Hv  += ash*Cv1
        Ha1 += ash*Ca1
        Ha2 += ash*Ca2
        Ha3 += ash*(Ca1+Ca3)
        Ht1 += ash*(Ct1+0.5*(w-1)*(Ct2-Ct3))
        Ht2 += ash*0.5*(w+1)*(Ct2+Ct3)
        Ht3 += ash*(Ct2)

        # Epsilon b and epsilon c
        eb = auxparams["eb"]
        ec = auxparams["ec"]

        IWs = self.IW(w, auxparams)
        Li1 = calc_Li1(w, IWs)
        
        # Hps += ec*(Li1[2]+Li1[3]*(w-1)+Li1[5]-Li1[6]*(w+1)) + eb*(Li1[1]-Li1[4])
        Hv  += ec*(Li1[2]-Li1[5]) + eb*(Li1[1]-Li1[4])
        Ha1 += ec*(Li1[2]-Li1[5]*(w-1)/(w+1)) + eb*(Li1[1]-Li1[4]*(w-1)/(w+1))
        Ha2 += ec*(Li1[3]+Li1[6])
        Ha3 += ec*(Li1[2]-Li1[3]-Li1[5]+Li1[6]) + eb*(Li1[1]-Li1[4])
        Ht1 += ec*(Li1[2]) + eb*(Li1[1])
        Ht2 += ec*(Li1[5]) - eb*(Li1[4])
        Ht3 += ec*(Li1[6]-Li1[3])

        # Upsilon Expansion
        corrb = eb*auxparams["upsilonb"]
        corrc = ec*auxparams["upsilonc"]
        # Hps += (corrc + corrb)
        Hv  += (corrc + corrb)
        Ha1 += (corrc + corrb)*(w-1)/(w+1)
        Ha2 += -2.*corrc/(w+1.)
        Ha3 += (corrc * (w-1.)/(w+1.)+corrb)
        Ht1 += 0.
        Ht2 += -(corrc-corrb)
        Ht3 += -2.*corrc/(w+1.)

        # Epsilon c squared
        ec2 = ec*ec
        Li2 = calc_Li2(w, IWs, auxparams["la2/laB2"])
        # Hps += ec2*(Li2[2]+Li2[3]*(w-1)+Li2[5]-Li2[6]*(w+1))
        Hv  += ec2*(Li2[2]-Li2[5])
        Ha1 += ec2*(Li2[2]-Li2[5]*(w-1)/(w+1))
        Ha2 += ec2*(Li2[3]+Li2[6])
        Ha3 += ec2*(Li2[2]-Li2[3]-Li2[5]+Li2[6])
        Ht1 += ec2*(Li2[2])
        Ht2 += ec2*(Li2[5])
        Ht3 += ec2*(Li2[6]-Li2[3])

        # Epsilon b squared
        if self.get_userparam("With1OverMb2"):
            eb2 = eb*eb
            # Hps += eb2*(Li2[1]-Li2[4])
            Hv  += eb2*(Li2[1]-Li2[4])
            Ha1 += eb2*(Li2[1]-Li2[4]*(w-1)/(w+1))
            Ha2 += 0.
            Ha3 += eb2*(Li2[1]-Li2[4])
            Ht1 += eb2*(Li2[1])
            Ht2 += -eb2*(Li2[4])
            Ht3 += 0.

        # Epsilon b * Epsilon c
        if self.get_userparam("With1OverMbMc"):
            eceb = ec*eb
            Mi = calc_Mi(w, IWs, auxparams["la1/laB2"], auxparams["la2/laB2"])
            Hv  += eceb * ((Mi[2]+Mi[9]) - (Mi[16]+Mi[18]))
            Ha1 += eceb * ((Mi[2]+Mi[9]) - (w-1.)/(w+1.)*(Mi[16]+Mi[18]))
            Ha2 += eceb * ((Mi[3]-Mi[10]) + (Mi[17]-Mi[19]))
            Ha3 += eceb * ((Mi[2]+Mi[9]) - (Mi[3]-Mi[10]) - (Mi[16]+Mi[18]) + (Mi[17]-Mi[19]))
            # Hps += eceb * (Mi[2] - Mi[9] + (w - 1.)*(Mi[3] + Mi[10]) + Mi[16] - Mi[18] - (1. + w )*(Mi[17]+Mi[19]))
            Ht1 += eceb * (Mi[2]-Mi[9])
            Ht2 += eceb * (Mi[16]-Mi[18])
            Ht3 += -eceb * ((Mi[3]+Mi[10]) - (Mi[17]+Mi[19]))

        # alpha_s/mb, alpha_s/mc
        if self.get_userparam("WithAsOverM"):
            asec = ash*ec
            aseb = ash*eb
            cmagc = auxparams["cmagc"]
            cmagb = auxparams["cmagb"]
            Cv2 = hqet.CV2(w, zBC)
            Cv3 = hqet.CV3(w, zBC)
            bound_lo = 1.0 + md._DER_EPS
            # dCps = md.derivative(lambda x: hqet.CP(x, zBC),  w, bound_lo = bound_lo)
            dCv1 = md.derivative(lambda x: hqet.CV1(x, zBC), w, bound_lo = bound_lo)
            dCa1 = md.derivative(lambda x: hqet.CA1(x, zBC), w, bound_lo = bound_lo)
            dCa2 = md.derivative(lambda x: hqet.CA2(x, zBC), w, bound_lo = bound_lo)
            dCa3 = md.derivative(lambda x: hqet.CA3(x, zBC), w, bound_lo = bound_lo)
            dCt1 = md.derivative(lambda x: hqet.CT1(x, zBC), w, bound_lo = bound_lo)
            dCt2 = md.derivative(lambda x: hqet.CT2(x, zBC), w, bound_lo = bound_lo)
            dCt3 = md.derivative(lambda x: hqet.CT3(x, zBC), w, bound_lo = bound_lo)

            Hv  += asec * (cmagc*Li1[2] + (Li1[2] - Li1[5])*Cv1 - (Li1[4] - Li1[5])*Cv3 + 2*(-1 + w)*dCv1)
            Ha1 += asec * (cmagc*Li1[2] + (Li1[2] - (Li1[5]*(-1 + w))/(1 + w))*Ca1 + ((Li1[4] - Li1[5])*(-1 + w)*Ca3)/(1 + w) + 2*(-1 + w)*dCa1)
            Ha2 += asec * (cmagc*Li1[3] + (Li1[3] + Li1[6])*Ca1 + (Li1[2] + Li1[5] + Li1[3]*(-1 + w) - Li1[6]*(1 + w))*Ca2 - ((Li1[4] - 3*Li1[5])*Ca3)/(1 + w) + 2*(-1 + w)*dCa2)
            Ha3 += asec * (cmagc*(-Li1[3] + Li1[2]) + (Li1[2] - Li1[3] - Li1[5] + Li1[6])*Ca1 + ((Li1[4] - 3*Li1[5])*w*Ca3)/(1 + w) + (Li1[2] + Li1[5] + Li1[3]*(-1 + w) - Li1[6]*(1 + w))*Ca3 + 2*(-1 + w)*(dCa1 + dCa3))
            # Hps += asec * ((cmagc*(Li1[3]*(-1 + w) + Li1[2])) + (Li1[2] + Li1[5] + Li1[3]*(-1 + w) - Li1[6]*(1 + w))*Cps + 2*(-1 + w)*dCps)
            Ht1 += asec * (cmagc*Li1[2] + Li1[2]*(Ct1 + ((-1 + w)*(Ct2 - Ct3))/2.) - (Li1[5]*(-1 + w)*(Ct2 - Ct3))/2. + Li1[5]*(-1 + w)*Ct3 + 2*(-1 + w)*(dCt1 + ((-1 + w)*(dCt2 - dCt3))/2.))
            Ht2 += asec * ((Li1[4] - Li1[5]*w)*Ct3 + (Li1[2]*(1 + w)*(Ct2 + Ct3))/2. + Li1[5]*(Ct1 - ((1 + w)*(Ct2 + Ct3))/2.) + (-1 + pow(w, 2))*(dCt2 + dCt3))
            Ht3 += asec * (-(cmagc*Li1[3]) - (Li1[3] - Li1[6])*Ct1 + (Li1[2] - Li1[5])*Ct2 - ((Li1[4] - 3*Li1[5])*Ct3)/(1 + w) + 2*(-1 + w)*dCt2)
            
            Hv  += aseb * (cmagb*Li1[1] + (Li1[1] - Li1[4])*Cv1 - (Li1[4] - Li1[5])*Cv2 + 2*(-1 + w)*dCv1)
            Ha1 += aseb * (cmagb*Li1[1] + (Li1[1] - (Li1[4]*(-1 + w))/(1 + w))*Ca1 + ((Li1[4] - Li1[5])*(-1 + w)*Ca2)/(1 + w) + 2*(-1 + w)*dCa1)
            Ha2 += aseb * ((Li1[1] - Li1[4] - 2*Li1[5])*Ca2 - ((Li1[4] - 3*Li1[5])*Ca2)/(1 + w) + 2*(-1 + w)*dCa2)
            Ha3 += aseb * (cmagb*Li1[1] + 2*Li1[5]*Ca2 + ((Li1[4] - 3*Li1[5])*w*Ca2)/(1 + w) + (Li1[1] - Li1[4])*(Ca1 + Ca3) + 2*(-1 + w)*(dCa1 + dCa3))
            # Hps += aseb * (cmagb*Li1[1] + (Li1[1] - Li1[4])*Cps + 2*(-1 + w)*dCps)
            Ht1 += aseb * (cmagb*Li1[1] - Li1[5]*(-1 + w)*Ct2 + Li1[1]*(Ct1 + ((-1 + w)*(Ct2 - Ct3))/2.) - (Li1[4]*(-1 + w)*(Ct2 - Ct3))/2. + 2*(-1 + w)*(dCt1 + ((-1 + w)*(dCt2 - dCt3))/2.))
            Ht2 += aseb * ((Li1[4] - Li1[5]*w)*Ct2 + (Li1[1]*(1 + w)*(Ct2 + Ct3))/2. - Li1[4]*(Ct1 + ((1 + w)*(Ct2 + Ct3))/2.) + (-1 + pow(w,2))*(dCt2 + dCt3))
            Ht3 += aseb * ((Li1[1] - Li1[4] - 2*Li1[5])*Ct2 - ((Li1[4] - 3*Li1[5])*Ct2)/(1 + w) + 2*(-1 + w)*dCt2)

        # alpha_s**2 for Ha1
        if self.get_userparam("WithHA1As2"):
            Ha1 += -0.944*(4./3.)*ash*ash

        h = {
            "V"  : Hv,
            "A1" : Ha1,
            "A2" : Ha2,
            "A3" : Ha3,
            "T1" : Ht1,
            "T2" : Ht2,
            "T3" : Ht3
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
            Mc = self.get_param(f"m_{self.V}")
            w = max(self.w(q2), 1)
            # IG
            Xi = self.xi(w)
            hhat = self.get_hhat(q2)
            h = {iff : Xi*hhat[iff] for iff in hhat}
            # NOTE: this performs the translation https://arxiv.org/pdf/1309.0301 eqns. 38-39,
            # should be analgous to eqns B7-B13 in https://arxiv.org/pdf/1908.09398
            return h_to_A(Mb, Mc, h, q2)