import flavio
from math import sqrt
from slophep.Predictions.FormFactorBase import FormFactor

class FormFactorBToV(FormFactor):
    def __init__(self, 
                 B: str,
                 V: str,
                 par: dict = None, 
                 scale: float = None):
        super().__init__(par, scale)
        
        self._B = B
        self._V = V
        self._name: str            = "BToV_FFBase"
        self._internalparams: dict = {
            "Mb"         : self.par[f'm_{self.B}'],
            "Mc"         : self.par[f'm_{self.V}'],
        }

    @property
    def B(self) -> str:
        """The B meson"""
        return self._B
    @property
    def V(self) -> str:
        """The Charmed Vector meson"""
        return self._V
    
    def w(self, q2: float) -> float:
        mB = self.internalparams["Mb"]
        mV = self.internalparams["Mc"]
        return (mB**2 + mV**2 - q2) / (2 * mB * mV)

    def get_ff(self, q2: float) -> dict:
        """Calculate form factors at particular q2. To implement in derived class.
        Must return in basis V, A0, A1, A12, T1, T2, T23

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            dictionary with FFs V, A0, A1, A12, T1, T2, T23
        """
        return {
            "V"   : 0.0,
            "A1"  : 0.0,
            "A0"  : 0.0,
            "T1"  : 0.0,
            "T2"  : 0.0,
            "A12" : 0.0,
            "T23" : 0.0
        }
    
    def get_ff_gfF1F2_basis(self, q2: float) -> dict:
        """Calculate FFs in common BGL basis g, f, F1, F2

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            dictionary with FFs g, f, F1, F2
        """
        Mb = self.internalparams["Mb"]
        Mc = self.internalparams["Mc"]

        ff_A = self.get_ff(q2)
        # Relations taken from https://github.com/eos/eos/blob/master/eos/form-factors/parametric-bgl1997-impl.hh lines 196 onwards
        ff = {
            "g"  : (2.0/(Mb + Mc))*ff_A["V"],
            "f"  : (Mb+Mc)*ff_A["A1"],
            "F1" : (8*Mb*Mc)*ff_A["A12"],
            "F2" : 2.0*ff_A["A0"]
        }
        return ff

    def get_ff_h_basis(self, q2: float) -> dict:
        """Calculate form factors in HQET basis h_i

        NOTE: hT3 and thus T3 must be computed from T23 and this requires
        dividing by a term that is zero at zero recoil. Because of this,
        the zero-recoil value is slightly offset from q2_max -> 0.999999999*q2_max
        
        This affects only this transformation of basis. The method
        get_ff(q2) does not have this issue.
        
        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            dictionary with FFs hV, hA0, hA1, hA2, hA3, hT1, hT2, hT3
        """
        
        mB = self.internalparams["Mb"]
        mV = self.internalparams["Mc"]
        # The kalman function lamdadstar is zero at zero recoil
        # To get T3 we need to divide by this and it breaks if lambdadstar is too small
        # This should only matter at the extreme edge of the q2 range and for this particular basis transformation    
        q2 = q2 if q2 < (mB - mV)**2 else 0.999999999*q2
        lambdabdstar = mB**4+mV**4+q2**2-2*(mB**2*mV**2+mB**2*q2+mV**2*q2)
        r = mV/mB
        sqrtr = sqrt(r)
        w = max((mB**2 + mV**2 - q2) / (2 * mB * mV), 1)

        ff_A = self.get_ff(q2)
        ff = {}
        ff["hV"] = 2*sqrtr/(1 + r) * ff_A["V"]
        ff["hA1"] = (1+r)/(sqrtr*(1+w)) * ff_A["A1"]
        A2 = -((16. * mB * mV**2 * (mB + mV))*ff_A["A12"] - ff_A["A1"] * (mB + mV)**2 * (mB**2 - mV**2 - q2))
        A2 /= (mB**4 + (mV**2 - q2)**2 - 2 * mB**2 * (mV**2 + q2))
        ff["hA2"] = (2*sqrtr*ff_A["A0"] - (r-w)*2*sqrtr/(1 + r)*A2
                      -(w+1)*ff["hA1"])/(2*r*w - r**2 - 1)
        ff["hA3"] = 2*sqrtr/(1 + r)*A2 - r*ff["hA2"]
        ff["hT1"] = (r+1)*((r-1)**2*ff_A["T2"] -2*r*(w-1)*ff_A["T1"]) / (2*sqrtr*(r**2 - 2*r*w + 1))
        ff["hT2"] = (1/(1-r))*((1+r)*ff["hT1"] - 2*sqrtr*ff_A["T1"])
        # Here is the problematic division by lambdastar at zero recoil
        T3 = 1/(lambdabdstar) * ((mB**2 - mV**2) * (mB**2 + 3 * mV**2 - q2) * ff_A['T2'] 
                                - (8 * mB * (mB - mV) * mV**2) * ff_A["T23"])
        ff["hT3"] = 1/(1-r**2)*(2*sqrtr*T3 - (1-r)*ff["hT1"] + (1+r)*ff["hT2"])

        return ff
    
    def get_ff_R_basis(self, q2: float) -> dict:
        """Calculate FF CLN ratios R0, R1, R2, formulae from ancillary material in https://arxiv.org/abs/2304.03137.
        Directly lifted from the ancillary files (LOAD_FIT.py)

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            dictionary with FFs hA1, R0, R1, R2
        """
        
        mB = self.internalparams["Mb"]
        mV = self.internalparams["Mc"]
        r = mV/mB
        ffh = self.get_ff_h(q2)
        hA1 = ffh["hA1"]
        hA2 = ffh["hA2"]
        hA3 = ffh["hA3"]
        hV = ffh["hV"]

        w = max((mB**2 + mV**2 - q2) / (2 * mB * mV), 1)
        R0=(1.0/(1.0+r))*(w +1 + w*(r*hA2-hA3)/hA1 - (hA2-r*hA3)/hA1)
        R1=hV/hA1
        R2=(r*hA2+hA3)/hA1
        return {"hA1" : hA1, "R0" : R0, "R1" : R1, "R2" : R2}