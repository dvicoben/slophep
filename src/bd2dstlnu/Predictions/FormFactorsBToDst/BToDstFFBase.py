import flavio
from math import sqrt

class FormFactor:
    def __init__(self, 
                 par: dict = None, 
                 scale: float = None):
        self._par: dict            = flavio.default_parameters.get_central_all() if type(par) == type(None) else par
        self._scale: float         = 4.8 if not scale else scale
        
        self._name: str            = "FFBase"
        self._ffpar: dict          = {}
        self._params: list         = []
        self._internalparams: dict = {
            "Mb"         : self.par['m_B0'],
            "Mc"         : self.par['m_D*+'],
        }

    @property
    def name(self) -> str: 
        """Name of FF scheme"""
        return self._name
    @property
    def scale(self) -> float:
        """Renorm scale"""
        return self._scale
    @property
    def par(self) -> dict:
        """Dictionary of flavio parameters, defaults to `flavio.default_parameters.get_central_all()`"""
        return self._par
    @property
    def ffpar(self) -> dict: 
        """Form factor (floating) parameters, make sure to overwrite in derived classes"""
        return self._ffpar
    @property
    def params(self) -> list[str]: 
        """List of (floating) form factor parameters (the ones in ffpar), make sure to overwrite in derived classes"""
        return self._params
    @property
    def internalparams(self) -> dict: 
        """Internal parameters used in FF computation (e.g. alpha_s, quark/meson masses), 
        where possible take from FormFactor.par, but things may differ if an implementation is taken from HAMMER
        so this is kept seperate"""
        return self._internalparams

    def _fullsetter(self, params: dict, constants: dict):
        """For usage with fluctuate - effectively alias of set_ff

        Parameters
        ----------
        params : dict
        """
        p = {**params, **constants}
        self.set_ff(**p)

    def set_ff(self, **kwargs):
        """Set form factors - argument names should be as in FormFactor.params"""
        for elem in kwargs:
            if elem not in self.ffpar:
                print(f"{elem} is not a parameter for FF {self.name}")
                continue
            if kwargs[elem] == None:
                continue
            self._ffpar[elem] = kwargs[elem]

    def set_ff_fromlist(self, params: list[float]):
        """Set form factors

        Parameters
        ----------
        params : list[float]
            Form factor parameters, in order as in FormFactor.params
        """
        pars = {self.params[k] : params[k] for k in range(len(self.params))}
        self.set_ff(**pars)

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
    
    def get_ff_gfF1F2(self, q2: float) -> dict:
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

    def get_ff_h(self, q2: float) -> dict:
        """Calculate form factors in HQET basis h_i

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
    
    def get_ff_R(self, q2: float) -> dict:
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