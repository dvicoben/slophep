import numpy as np
from slophep.Predictions.FormFactorsBToP import FormFactorBToP
from flavio.physics.bdecays.formfactors import hqet
from flavio.physics.bdecays.formfactors.b_p.cln import h_to_f
import slophep.Predictions.Math.derivatives as md


class BLPRXP_BToP(FormFactorBToP):
    def __init__(self, B: str, P: str, par: dict = None, scale: float = None, *ffargs):
        super().__init__(B, P, par, scale)
        
        self._name = "BToP_BLPRXP"
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
            "WithAsOverM"   : True
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
        return (np.sqrt(w+1) - sqrt2*a)/(np.sqrt(w+1) + sqrt2*a)

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
    
    def IW(self, w: float) -> dict:
        IWs = {
            "CHI2"  : self.ffpar["Chi21"] + (w-1.0)*self.ffpar["Chi2p"],
            "CHI3"  : self.ffpar["Chi3p"]*(w-1.0),
            "ETA"   : self.ffpar["Eta1"] + (w-1.0)*self.ffpar["Etap"],
            "BETA1" : self.internalparams["la1/laB2"]/4.0,
            "BETA2" : self.ffpar["Beta21"],
            "BETA3" : self.internalparams["la2/laB2"]/8.0 + self.ffpar["Beta3p"]*(w-1.0),
            "PHI1"  : (self.internalparams["la1/laB2"]/3. - self.internalparams["la2/laB2"]/2.)/2. + (self.ffpar["Phi1p"])*(w-1.),
            "PHI1Q" : self.ffpar["Phi1p"]
        }
        return IWs

    def Li1(self, w: float, IWs: dict = None) -> list[float]:
        IWs = IWs if IWs is not None else self.IW(w)
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
    
    def Li2(self, w: float, IWs: dict = None) -> list[float]:
        IWs = IWs if IWs is not None else self.IW(w)
        la2OverlaB2 = self.internalparams["la2/laB2"]
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

    def Mi(self, w: float, IWs: dict = None) -> list[float]:
        IWs = IWs if IWs is not None else self.IW(w)
        la1OverlaB2 = self.internalparams["la1/laB2"]
        la2OverlaB2 = self.internalparams["la2/laB2"]
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

    def get_hhat(self, q2: float) -> dict:
        w = max(self.w(q2), 1)
        zBC = self.internalparams["zBC"]

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
        ash = self.internalparams["ash"]
        Hp += ash*(Cv1+0.5*(w+1)*(Cv2+Cv3))
        Hm += ash*0.5*(w+1)*(Cv2-Cv3)
        Ht += ash*(Ct1-Ct2+Ct3)

        # Epsilon b and epsilon c
        eb = self.internalparams["eb"]
        ec = self.internalparams["ec"]

        IWs = self.IW(w)
        Li1 = self.Li1(w, IWs)
        
        Hp += (ec + eb)*Li1[1]
        Hm += (ec - eb)*Li1[4]
        Ht += (ec + eb)*(Li1[1]-Li1[4])

        # Upsilon Expansion
        corrb = eb*self.internalparams["upsilonb"]
        corrc = ec*self.internalparams["upsilonc"]
        Hp += 0.
        Hm += -(corrc-corrb)
        Ht += (corrc+corrb)

        # Epsilon c squared
        ec2 = ec*ec
        Li2 = self.Li2(w, IWs)
        # Hps += ec2*(Li2[2]+Li2[3]*(w-1)+Li2[5]-Li2[6]*(w+1))
        Hp += ec2*Li2[1]
        Hm += ec2*Li2[4]
        Ht += ec2*(Li2[1]-Li2[4])

        # Epsilon b squared
        if self.internalparams["With1OverMb2"]:
            eb2 = eb*eb
            Hp += +eb2*Li2[1]
            Hm += -eb2*Li2[4]
            Ht += +eb2*(Li2[1]-Li2[4])

        # Epsilon b * Epsilon c
        if self.internalparams["With1OverMbMc"]:
            eceb = ec*eb
            Mi = self.Mi(w, IWs)
            Hp += eceb * (Mi[1] - Mi[8])
            Hm += 0.
            Ht += eceb * (Mi[1] + Mi[8] - 2.*Mi[15])

        # alpha_s/mb, alpha_s/mc
        if self.internalparams["WithAsOverM"]:
            asec = ash*ec
            aseb = ash*eb
            cmagc = self.internalparams["cmagc"]
            cmagb = self.internalparams["cmagb"]
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
        # IG
        Xi = self.xi(w)
        hhat = self.get_hhat(q2)
        h = {iff : Xi*hhat[iff] for iff in hhat}
        
        return h_to_f(Mb, Mc, h, q2)
