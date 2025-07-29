from math import sqrt
import numpy as np
from slophep.Predictions.FormFactorsBToV import FormFactorBToV


class BSZ_BToV(FormFactorBToV):
    def __init__(self, B: str, V: str, par: dict = None, scale: float = None, *ffargs):
        super().__init__(B, V, par, scale)
        
        self._name = "BSZ"
        # Parameters from EOS, see https://eoshep.org/doc/reference/parameters.html#parameters-in-b-to-v-form-factor-parametrizations
        # Should correspond to LCSR+Lattice fit from arXiv:1811.00983v2
        self._ffpar = { 
            "A0_0"   : 0.6759274184148141,
            "A0_1"   : -2.6637000799846655,
            "A0_2"   : 9.24495323738047,
            "A1_0"   : 0.601315572767782,
            "A1_1"   : 0.2689438444495897,
            "A1_2"   : 10.486032487112988,
            "A12_1"  : -0.14876383294863582,
            "A12_2"  : 5.724359813703122,
            "V_0"    : 0.6851323923652662,
            "V_1"    : -1.1678379355866901,
            "V_2"    : -4.016847943717434,
            "T1_0"   : 0.6302298493900512,
            "T1_1"   : -1.4175767766320952,
            "T1_2"   : -0.6878528871949201,
            "T2_1"   : 1.5226848746033717,
            "T2_2"   : 0.006625613047247359,
            "T23_0"  : 0.806230464649646,
            "T23_1"  : 0.5855443911134184,
            "T23_2"  : 4.667079193615467,
        }
        self._params = ["A0_0", "A0_1", "A0_2", 
                        "A1_0", "A1_1", "A1_2", 
                        "A12_1", "A12_2",
                        "V_0", "V_1", "V_2",
                        "T1_0", "T1_1", "T1_2",
                        "T2_1", "T2_2",
                        "T23_0", "T23_1", "T23_2"]
        internalparams = {
            "m0"  : 6.275,
            "m1m" : 6.330,
            "m1p" : 6.767
        }
        self._internalparams.update(internalparams)

    def z(self, mB: float, mM: float, q2: float, t0: float = None):
        """Form factor expansion parameter z.

        Parameters
        ----------

        mB : float
            initial pseudoscalar meson mass
        mM : float
            final meson meson mass
        q2 : float
            momentum transfer squared q2
        t0 : float
            parameter t_0.
            If not given, chosen as t_0 = t_+ (1 - sqrt(1 - t\_-/t_+))$ where
            t_+- = (m_B +- m_M)^2.
            If equal to `tm`, set to t_0=t\_-
        """
        tm = (mB-mM)**2
        tp = (mB+mM)**2
        if t0 is None:
            t0 = tp*(1-sqrt(1-tm/tp))
        elif t0 == 'tm':
            t0 = tm
        sq2 = sqrt(tp-q2)
        st0 = sqrt(tp-t0)
        return (sq2-st0)/(sq2+st0)
    
    def zs(self, mB: float, mV: float, q2: float, t0: float = None):
        zq2 = self.z(mB, mV, q2, t0)
        z0 = self.z(mB, mV, 0, t0)
        return np.array([1, zq2-z0, (zq2-z0)**2])
    
    def pole(self, q2: float, mres: float):
        return 1/(1-q2/(mres**2))

    def get_ff(self, q2: float) -> dict:
        """Calculates BSZ form factors, following https://arxiv.org/pdf/1503.05534 and https://arxiv.org/pdf/1811.00983.
        Implementation lifted from flavio https://github.com/flav-io/flavio/blob/master/flavio/physics/bdecays/formfactors/b_v/bsz.py, 
        and should be in line with EOS https://github.com/eos/eos/blob/master/eos/form-factors/parametric-bsz2015-impl.hh

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            FF dictionary
        """
        # Resonance masses used in arXiv:1811.00983
        m0  = self.internalparams["m0"]
        m1m = self.internalparams["m1m"]
        m1p = self.internalparams["m1p"]
        # mres = [m0, m1m, m1p]

        mB = self.internalparams["Mb"]
        mV = self.internalparams["Mc"]


        aA0 = np.array([self.ffpar["A0_0"], self.ffpar["A0_1"], self.ffpar["A0_2"]])
        aA1 = np.array([self.ffpar["A1_0"], self.ffpar["A1_1"], self.ffpar["A1_2"]])
        aA12_0 = self.ffpar["A0_0"] / (8*mB*mV / (mB**2-mV**2)) # https://arxiv.org/pdf/1503.05534 eqn. 17
        aA12 = np.array([aA12_0, self.ffpar["A12_1"], self.ffpar["A12_2"]])
        aV = np.array([self.ffpar["V_0"], self.ffpar["V_1"], self.ffpar["V_2"]])
        aT1 = np.array([self.ffpar["T1_0"], self.ffpar["T1_1"], self.ffpar["T1_2"]])
        aT2 = np.array([self.ffpar["T1_0"], self.ffpar["T2_1"], self.ffpar["T2_2"]]) # https://arxiv.org/pdf/1503.05534 eqn. 17 T2_0 = T1_0
        aT23 = np.array([self.ffpar["T23_0"], self.ffpar["T23_1"], self.ffpar["T23_2"]])

        # mresdict = {'A0': 0,'A1': 2,'A12': 2,'V': 1,'T1': 1,'T2': 2,'T23': 2}
        z_exp = self.zs(mB, mV, q2, t0=None)
        ff = {
            "A0"  : self.pole(q2, m0)*np.dot(aA0, z_exp),
            "A1"  : self.pole(q2, m1p)*np.dot(aA1, z_exp),
            "A12" : self.pole(q2, m1p)*np.dot(aA12, z_exp),
            "V"   : self.pole(q2, m1m)*np.dot(aV, z_exp),
            "T1"  : self.pole(q2, m1m)*np.dot(aT1, z_exp),
            "T2"  : self.pole(q2, m1p)*np.dot(aT2, z_exp),
            "T23" : self.pole(q2, m1p)*np.dot(aT23, z_exp),
        }
        
        return ff
