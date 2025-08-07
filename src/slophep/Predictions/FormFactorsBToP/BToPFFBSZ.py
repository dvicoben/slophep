from math import sqrt
import numpy as np
from slophep.Predictions.FormFactorsBToP import FormFactorBToP

class BSZ_BToP(FormFactorBToP):
    def __init__(self, B: str, P: str, par: dict = None, scale: float = None, *ffargs):
        super().__init__(B, P, par, scale)

        self._name = "BSZ"
        self._ffpar = {

        }
        # Parameters from B->P EOS https://eoshep.org/doc/reference/parameters.html#parameters-in-b-to-p-form-factor-parametrizations
        self._ffpar = {
            "f+_0" : 0.680881482963591,
            "f+_1" : -2.221782088112512,
            "f+_2" : 0.6541000835887975,
            "f0_1" : 0.5779088408438552,
            "f0_2" : 1.3271162015486215,
            "fT_0" : 0.5648201165054019,
            "fT_1" : -2.5027753865334676,
            "fT_2" : 10.591752810763625
        }
        self._params = [
            "f+_0", "f+_1", "f+_2",
            "f0_1", "f0_2",
            "fT_0", "fT_1", "fT_2"
        ]
        internalparams = {
            "m0" : 6.420,
            "mp" : 6.330
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
    
    def zs(self, mB: float, mM: float, q2: float, t0: float = None):
        zq2 = self.z(mB, mM, q2, t0)
        z0 = self.z(mB, mM, 0, t0)
        return np.array([1, zq2-z0, (zq2-z0)**2])
    
    def pole(self, q2: float, mres: float):
        return 1/(1-q2/(mres**2))
    
    def get_ff(self, q2: float) -> dict:
        """Calculates BSZ form factors, following https://arxiv.org/pdf/1503.05534 and https://arxiv.org/pdf/1811.00983.
        Implementation lifted from flavio https://github.com/flav-io/flavio/blob/master/flavio/physics/bdecays/formfactors/b_p/bsz.py, 

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            FF dictionary
        """
        # Resonance masses used in arXiv:1811.00983
        mpl = self.internalparams["mp"]
        m0 = self.internalparams["m0"]
        mB = self.internalparams["Mb"]
        mP = self.internalparams["Mc"]

        afp = np.array([self.ffpar["f+_0"], self.ffpar["f+_1"], self.ffpar["f+_2"]])
        # f0_0 is fixed by kinematic constraint f0(0) = f+(0) so f0_0 = f+_0
        af0 = np.array([self.ffpar["f+_0"], self.ffpar["f0_1"], self.ffpar["f0_2"]])
        afT = np.array([self.ffpar["fT_0"], self.ffpar["fT_1"], self.ffpar["fT_2"]])

        z_exp = self.zs(mB, mP, q2, t0=None)
        ff = {
            "f+" : self.pole(q2, mpl)*np.dot(afp, z_exp),
            "f0" : self.pole(q2, m0)*np.dot(af0, z_exp),
            "fT" : self.pole(q2, mpl)*np.dot(afT, z_exp)
        }
        
        return ff