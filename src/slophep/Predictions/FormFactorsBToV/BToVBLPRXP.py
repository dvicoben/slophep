from math import sqrt
import numpy as np
from slophep.Predictions.FormFactorsBToV import FormFactorBToV
from flavio.physics.bdecays.formfactors import hqet
from flavio.physics.bdecays.formfactors.b_v.cln import h_to_A
from slophep.Predictions.Math.derivatives import derivative


class BLPRXP_BToV(FormFactorBToV):
    def __init__(self, B: str, V: str, par: dict = None, scale: float = None, *ffargs):
        super().__init__(B, V, par, scale)
        
        self._name = "BToV_BLPRXP"
        self._ffpar = {
            "RhoStSq" : 1.104,
            "cSt"     : 2.392,
            "Chi21"   : -0.116,
            "Chi2p"   : 0.0,
            "Chi3p"   : 0.0,
            "Eta1"    : 0.336,
            "Etap"    : 0.0,
            "Phi1p"   : 0.252,
            "Beta21"  : 0.0,
            "Beta3p"  : 0.0
        }
        # "delta_RhoSq","delta_cSt","delta_chi21","delta_chi2p","delta_chi3p","delta_eta1","delta_etap","delta_phi1p","delta_beta21","delta_beta3p"
        self._params = ["RhoStSq", "cSt",
                        "Chi21", "Chi2p", "Chi3p", "Eta1", "Etap", 
                        "Phi1p", "Beta21", "Beta3p"]
        
        self.calculate_internal()

    def calculate_internal(self, rawinternalpar: dict = {}):
        # 1S scheme: mb1S = 4710, mbBar = 5313, lambda1 = -3 10^5 MeV^2, delta mb-mc = 3400, alpha_s = 26/100
        params = {
            "a"             : 1.509/np.sqrt(2.0),
            "as"            : 0.27,
            "mb"            : 4.707,
            "mc"            : 4.707 - 3.407,
            "mBBar"         : 5.313,
            "mDBar"         : 1.973,
            "rho1"          : -0.363,
            "la2"           : 0.12,
            "rD"            : self.par['m_D0']/self.par['m_B0'],
            "With1OverMb2"  : True,
            "With1OverMbMc" : True,
            "WithAsOverM"   : True,
            "WithHA1As2"    : True
        }
        params.update(rawinternalpar)

        _as    = params["as"]
        _ash   = _as/np.pi
        _mb    = params["mb"]
        _mc    = params["mc"]
        _mBBar = params["mBBar"]
        _mDBar = params["mDBar"]
        _rho1  = params["rho1"]
        _la2   = params["la2"]
        # Some additional calculations to avoid recomputing:        
        zBC = _mc/_mc
        laB = (_mb*_mBBar - _mc*_mDBar)/(_mb-_mc) - (_mb+_mc) + _rho1/(4.*_mb*_mc)
        la1 = (2.*_mb*_mc)/(_mb - _mc)*(_mBBar - _mDBar - (_mb-_mc)) + _rho1*(_mb+_mc)/(2.*_mb*_mc)
        # _eb = laB/(2.*_mb)
        # _ec = laB/(2.*_mc)
        # _la2OverlaB2 = _la2/(laB*laB)
        # _la1OverlaB2 = la1/(laB*laB)
        corr1S = 2.* pow(_as/3. * 0.796, 2.)
        dmbc = _mb - _mc
        mb_1 = _mb * corr1S
        mc_1 = _mb * corr1S
        LambdaBar_0 = (_mb*_mBBar - _mc*_mDBar)/dmbc - (2.*_mb-dmbc) + _rho1/(4*_mb*_mc)
        LambdaBar_1 = (mb_1 * _mBBar - mc_1 * _mDBar)/dmbc - 2.*mb_1
        _upsilonb = LambdaBar_1/LambdaBar_0 - mb_1/_mb
        _upsilonc = LambdaBar_1/LambdaBar_0 - mc_1/_mc

        internalparams = {
            "ash"      : _ash,
            "zBC"      : zBC,
            "cmagb"    : 3./2.*(0.5*np.log(zBC) + 13./9.),
            "cmagc"    : -3./2.*(0.5*np.log(zBC) + 13./9.),
            "laB"      : laB,
            "la1"      : la1,
            "la2/laB2" : _la2/(laB*laB),
            "la1/laB2" : la1/(laB*laB),
            "corr1S"   : corr1S,
            "dmbc"     : dmbc,
            "mb1"      : mb_1,
            "mc1"      : mc_1,
            "LamBar0"  : LambdaBar_0,
            "LamBar1"  : LambdaBar_1,
            "eb"       : laB/(2.*_mb),
            "ec"       : laB/(2.*_mc),
            "upsilonb" : _upsilonb,
            "upsilonc" : _upsilonc
        }
        self.internalparams.update(params)
        self.internalparams.update(internalparams)

    def z_of_w(self, w: float, a: float) -> float:
        sqrt2 = np.sqrt(2)
        return (np.sqrt(w+1) - sqrt2*a)/(sqrt(w+1) + sqrt2*a)

    def xi(self, w: float) -> float:
        a = self.internalparams["a"]
        a2 = a**2
        a4 = a**4
        zCon = self.z_of_w(w, a)
        zCon1 = (1.-a)/(1.+a)
        RhoSq = self.ffpar["RhoStSq"]
        cSt = self.ffpar["cSt"]
        Xi = (
            (1. - 8*a2*RhoSq*zCon + 16.*(2*cSt * a4 - RhoSq * a2)*zCon*zCon)
            /(1. - 8*a2*RhoSq*zCon1 + 16.*(2*cSt * a4 - RhoSq * a2)*zCon1*zCon1)
        )
        return Xi
    
    def Li1(self, w: float, IWs: dict[str, float]) -> list[float]:
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

    def get_ff(self, q2: float) -> dict:
        """FF in BLPRXP parameterisation from ARXIV REF as in HAMMER v1.4.1

        Parameters
        ----------
        q2 : float
            q2 value to calculate FF at

        Returns
        -------
        dict
            FF dictionary
        """
        Mb = self.internalparams["Mb"]
        Mc = self.internalparams["Mc"]

        w = max(self.w(q2), 1)
        zBC = self.internalparams["zBC"]

        # a = self.internalparams["a"]
        # a2 = a**2
        # a4 = a**4
        # zCon = self.z_of_w(w, a)
        # zCon1 = (1.-a)/(1.+a)
        # RhoSq = self.ffpar["RhoStSq"]
        # cSt = self.ffpar["cSt"]
        # Xi = (
        #     (1. - 8*a2*RhoSq*zCon + 16.*(2*cSt * a4 - RhoSq * a2)*zCon*zCon)
        #     /(1. - 8*a2*RhoSq*zCon1 + 16.*(2*cSt * a4 - RhoSq * a2)*zCon1*zCon1)
        # )
        Xi = self.xi(w)

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
        ash = self.internalparams["ash"]
        # Hps += ash*Cps
        Hv  += ash*Cv1
        Ha1 += ash*Ca1
        Ha2 += ash*Ca2
        Ha3 += ash*(Ca1+Ca3)
        Ht1 += ash*(Ct1+0.5*(w-1)*(Ct2-Ct3))
        Ht2 += ash*0.5*(w+1)*(Ct2+Ct3)
        Ht3 += ash*(Ct2)


        # Epsilon b and epsilon c
        eb = self.internalparams["eb"]
        ec = self.internalparams["ec"]

        IWs = {
            "CHI2"  : self.ffpar["Chi21"] + (w-1.0)*self.ffpar["Chi2p"],
            "CHI3"  : self.ffpar["Chi3"]*(w-1.0),
            "ETA"   : self.ffpar["Eta1"] + (w-1.0)*self.ffpar["Etap"],
            "BETA1" : self.internalparams["la1/laB2"]/4.0,
            "BETA2" : self.ffpar["Beta21"],
            "BETA3" : self.internalparams["la2/laB2"]/8.0 + self.ffpar["Beta3p"]*(w-1.0),
            "PHI1"  : (self.internalparams["la1/laB2"]/3. - self.internalparams["la2/laB2"]/2.)/2. + (self.ffpar["Phi1p"])*(w-1.),
            "PHI1Q" : self.ffpar["Phi1p"]
        }
        IWsCHI2  = self.ffpar["Chi21"] + (w-1.0)*self.ffpar["Chi2p"]
        IWsCHI3  = self.ffpar["Chi3"]*(w-1.0)
        IWsETA   = self.ffpar["Eta1"] + (w-1.0)*self.ffpar["Etap"]
        IWsBETA1 = self.internalparams["la1/laB2"]/4.0
        IWsBETA2 = self.ffpar["Beta21"]
        IWsBETA3 = self.internalparams["la2/laB2"]/8.0 + self.ffpar["Beta3p"]*(w-1.0)
        IWsPHI1  = (self.internalparams["la1/laB2"]/3. - self.internalparams["la2/laB2"]/2.)/2. + (self.ffpar["Phi1p"])*(w-1.)
        IWsPHI1Q = self.ffpar["Phi1p"]
        
        Li1 = [
            0.0,
            4.*(3.* IWsCHI3 - (w-1.)*IWsCHI2),
            -4.*IWsCHI3,
            4.*IWsCHI2,
            2.* IWsETA - 1.,
            -1.0,
            -2.* (IWsETA + 1.)/(w + 1.)
        ]
        Hps += ec*(Li1[2]+Li1[3]*(w-1)+Li1[5]-Li1[6]*(w+1)) + eb*(Li1[1]-Li1[4])
        Hv  += ec*(Li1[2]-Li1[5]) + eb*(Li1[1]-Li1[4])
        Ha1 += ec*(Li1[2]-Li1[5]*(w-1)/(w+1)) + eb*(Li1[1]-Li1[4]*(w-1)/(w+1))
        Ha2 += ec*(Li1[3]+Li1[6])
        Ha3 += ec*(Li1[2]-Li1[3]-Li1[5]+Li1[6]) + eb*(Li1[1]-Li1[4])
        Ht1 += ec*(Li1[2]) + eb*(Li1[1])
        Ht2 += ec*(Li1[5]) - eb*(Li1[4])
        Ht3 += ec*(Li1[6]-Li1[3])

        # Upsilon Expansion
        corrb = eb*self.internalparams["upsilonb"]
        corrc = ec*self.internalparams["upsilonc"]
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
        la1OverlaB2 = self.internalparams["la1/laB2"]
        la2OverlaB2 = self.internalparams["la2/laB2"]
        Li2 = [
            0.0,
            2.*IWsBETA1 + 4.*(3.*IWsBETA3 - (w-1.)*IWsBETA2),
            2.*IWsBETA1 - 4.*IWsBETA3,
            4.*IWsBETA2,
            3*la2OverlaB2 + 2.*(w+1.)*IWsPHI1,
            la2OverlaB2 + 2.*(w+1.)*IWsPHI1,
            4.*IWsPHI1,
        ]
        # Hps += ec2*(Li2[2]+Li2[3]*(w-1)+Li2[5]-Li2[6]*(w+1))
        Hv  += ec2*(Li2[2]-Li2[5])
        Ha1 += ec2*(Li2[2]-Li2[5]*(w-1)/(w+1))
        Ha2 += ec2*(Li2[3]+Li2[6])
        Ha3 += ec2*(Li2[2]-Li2[3]-Li2[5]+Li2[6])
        Ht1 += ec2*(Li2[2])
        Ht2 += ec2*(Li2[5])
        Ht3 += ec2*(Li2[6]-Li2[3])

        # Epsilon b squared
        if self.internalparams["With1OverMb2"]:
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
        if self.internalparams["With1OverMbMc"]:
            eceb = ec*eb
            Mi = np.zeros(25) 
            Mi[8] = (
                la1OverlaB2+ 6.* la2OverlaB2/(w + 1.)
                -2.*(w-1.)*IWsPHI1 
                - 2.*(2.*IWsETA - 1.)*(w - 1.)/(w + 1.)
            )
            Mi[9] = 3.*la2OverlaB2/(w+1.) + 2.*IWsPHI1 - (2.*IWsETA - 1.) * (w - 1.)/(w + 1.)
            Mi[10] = (
                la1OverlaB2/3. 
                - la2OverlaB2*(w+4.)/(2.*(w+1.))
                + 2.*(w + 2.)*IWsPHI1Q - (2.*IWsETA - 1.)/(w + 1.)
            )
            Hv  += eceb * ((Mi[2]+Mi[9]) - (Mi[16]+Mi[18]))
            Ha1 += eceb * ((Mi[2]+Mi[9]) - (w-1.)/(w+1.)*(Mi[16]+Mi[18]))
            Ha2 += eceb * ((Mi[3]-Mi[10]) + (Mi[17]-Mi[19]))
            Ha3 += eceb * ((Mi[2]+Mi[9]) - (Mi[3]-Mi[10]) - (Mi[16]+Mi[18]) + (Mi[17]-Mi[19]))
            # Hps += eceb * (Mi[2] - Mi[9] + (w - 1.)*(Mi[3] + Mi[10]) + Mi[16] - Mi[18] - (1. + w )*(Mi[17]+Mi[19]))
            Ht1 += eceb * (Mi[2]-Mi[9])
            Ht2 += eceb * (Mi[16]-Mi[18])
            Ht3 += -eceb * ((Mi[3]+Mi[10]) - (Mi[17]+Mi[19]))

        # alpha_s/mb, alpha_s/mc
        if self.internalparams["WithAsOverM"]:
            asec = ash*ec
            aseb = ash*eb
            cmagc = self.internalparams["cmagc"]
            cmagb = self.internalparams["cmagb"]
            Cv2 = hqet.CV2(w, zBC)
            Cv3 = hqet.CV3(w, zBC)
            # dCps = derivative(lambda x: hqet.CP(x, zBC),  w, bound_lo = 1.0)
            dCv1 = derivative(lambda x: hqet.CV1(x, zBC), w, bound_lo = 1.0)
            dCa1 = derivative(lambda x: hqet.CA1(x, zBC), w, bound_lo = 1.0)
            dCa2 = derivative(lambda x: hqet.CA2(x, zBC), w, bound_lo = 1.0)
            dCa3 = derivative(lambda x: hqet.CA3(x, zBC), w, bound_lo = 1.0)
            dCt1 = derivative(lambda x: hqet.CT1(x, zBC), w, bound_lo = 1.0)
            dCt2 = derivative(lambda x: hqet.CT2(x, zBC), w, bound_lo = 1.0)
            dCt3 = derivative(lambda x: hqet.CT3(x, zBC), w, bound_lo = 1.0)

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


        if self.internalparams["WithHA1As2"]:
            Ha1 += -0.944*(4./3.)*ash*ash
        


        h = {
            "V"  : Xi*Hv,
            "A1" : Xi*Ha1,
            "A2" : Xi*Ha2,
            "A3" : Xi*Ha3,
            "T1" : Xi*Ht1,
            "T2" : Xi*Ht2,
            "T3" : Xi*Ht3
        }
        # NOTE: this performs the translation https://arxiv.org/pdf/1309.0301 eqns. 38-39,
        # should be analgous to eqns B7-B13 in https://arxiv.org/pdf/1908.09398
        return h_to_A(Mb, Mc, h, q2)

    def get_hhat(self, q2: float) -> dict:
        
        hhat = {
            # "V"  : Hv,
            # "A1" : Ha1,
            # "A2" : Ha2,
            # "A3" : Ha3,
            # "T1" : Ht1,
            # "T2" : Ht2,
            # "T3" : Ht3
        }
        return hhat